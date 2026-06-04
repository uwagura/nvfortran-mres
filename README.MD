# nvfortran Minimally Reproducible Examples (MREs)

A collection of MREs for nvfortran compiler bugs and unexpected behaviours found while
porting MOM6 to GPU with NVIDIA HPC SDK (nvfortran). Each MRE isolates a single issue
to the smallest possible code that reproduces it.

## Environment

All MREs use nvfortran 26.3. Activate the compiler before building:

```bash
source envs/linux-nvidia_2026_3.env
```

Build with `make <target>` from this directory. See `Makefile` for all targets.

---

## Issue 1 — Derived types failing to inline (nvfortran < 26.3)

**File:** `MRE_derived_types_failing_to_inline.F90`

Subroutines with derived-type arguments failed to inline at all in older nvfortran versions,
with no explanation. Fixed in nvfortran 26.3. The MRE demonstrates this by compiling with
`-DUSE_DERIVED_TYPES` and `-Minline=reshape -Minfo=inline`:

```bash
source envs/linux-nvidia_2026_3.env
nvfortran -O2 -mp=gpu -stdpar=gpu -gpu=mem:separate -Mnofma -r8 \
    -Minfo=inline -Minline=reshape -DUSE_DERIVED_TYPES \
    derived_types_mod.F90 MRE_derived_types_failing_to_inline.F90 -o mre_dt.exe
```

---

## Issue 2 — gradKE: "array reshaping not enabled" inlining failure

**File:** `MRE_gradKE_failing_inline_MOM6.F90`  
**Makefile targets:** `mre_gradke`, `mre_gradke_reshape`

A simplified analog of `gradKE` from `MOM_CoriolisAdv.F90`. When called with a 3-D array
slice (e.g., `u(:,:,k)`) as the argument for a 2-D dummy array, nvfortran requires the
`reshape` sub-option to `-Minline`. Without it, the compiler prints:

```
subprogram not inlined -- array reshaping not enabled: gradke_mre2, argument 1
```

**Reproduce:**
```bash
make mre_gradke          # shows "array reshaping not enabled"
make mre_gradke_reshape  # shows "gradke_mre2 inlined, size=..."
```

---

## Issue 3 — gradKE: silent non-inlining when derived type has a pointer to an external module type

**Files:** `external_type_mod.F90`, `MRE_external_ptr_silent_noinline.F90`  
**Makefile targets:** `mre_extptr_reshape`, `mre_extptr`

### Background

In the real MOM6 build, `gradKE` in `MOM_CoriolisAdv.F90` is still **not inlined** even
when both `-Minline=name:gradke,reshape` flags are passed, and the compiler produces
**absolutely no diagnostic output** at the call site — no "inlined" message, no
"subprogram not inlined" message. The surrounding `do concurrent` blocks do generate GPU
output, so the silence is specific to the `gradKE` call.

### Root cause

The trigger is: **a derived-type argument whose type definition contains a `POINTER` member
pointing to a type defined in a separately compiled module**.

In MOM6, `ocean_grid_type` (the type of the `G` argument to `gradKE`) contains:
```fortran
type(MOM_domain_type), pointer :: Domain => NULL()
```
where `MOM_domain_type` is defined in `MOM_domains` — a separately compiled module.
This single pointer member causes nvfortran's inliner to silently abandon the inline
attempt with zero diagnostic output.

The content and complexity of the pointed-to type are irrelevant — a pointer to a type
with just three integers in a separate module produces the same silence.

### What was ruled out

The following features, tested individually, do **not** cause the silent failure:
- Pointer members to locally-defined types (same compilation unit)
- Large numbers (40+) of allocatable array members
- Embedded non-pointer derived types from other modules
- Character members, 1D pointer arrays
- Allocatable arrays with fixed lower bounds (`dimension(0:,:)`)
- Self-referential pointer types (defined locally)

### The MRE

`external_type_mod.F90` defines a trivial 3-integer type in a separate module.
`MRE_external_ptr_silent_noinline.F90` contains two side-by-side variants:

| Variant | Grid type | Expected compiler output |
|---------|-----------|--------------------------|
| A — `grid_no_ptr` | No pointer to external type | `gradke_no_ptr inlined, size=45, ...` |
| B — `grid_with_ptr` | Has `type(external_type), pointer :: Domain` | *(complete silence)* |

**Reproduce:**
```bash
make mre_extptr_reshape
```

Expected `-Minfo=inline` output:
```
gradke_no_ptr:
caller_no_ptr:
     91, gradke_no_ptr inlined, size=45, file MRE_external_ptr_silent_noinline.F90 (63)
gradke_with_ptr:
caller_with_ptr:
                    ← no message here at all
```

---

## Other MREs

| File | Issue |
|------|-------|
| `MRE_subroutine_call_in_teams_region.F90` | Subroutine call inside `!$omp target teams` region |
| `mre1_llc_undefined_value.F90` | LLC "use of undefined value '%g'" compile error with unmapped derived type in `do concurrent` |
| `mre1a_llc_no_omp_parallel.F90` | Variant of MRE1 without `!$omp parallel do` |
| `mre2_do_concurrent_local.F90` | `do concurrent local()` segfault when variable also in `target enter data` map |
| `MRE_pointer_member_silent_noinline.F90` | Earlier investigation: 4-variant test of pointer members (all locally defined); all inline — showing the issue requires an *external* module type |
| `MRE_2d_array_privatization.F90` | 2-D array privatization in GPU `do concurrent` |
