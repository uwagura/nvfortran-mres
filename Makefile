FC = nvfortran
FFLAGS = -O2 -mp=gpu -stdpar=gpu -gpu=mem:separate -Mnofma -Minline=reshape -DUSE_DERIVED_TYPES
#FFLAGS = -O2 -mp=gpu -stdpar=gpu -gpu=mem:separate -Mnofma -Minfo=all -Minline=reshape

SRC = derived_types_mod.F90 MRE_subroutine_call_in_teams_region.F90 chain_mod4.F90 chain_mod3.F90 chain_mod2.F90 chain_mod1.F90 MRE_subroutine_call_in_teams_region_program.F90
OBJ = $(SRC:.F90=.o)
EXE = many_calls_with_derived_types_O2.exe

# MRE flags matching the MOM6 GPU build flags from ocean_only/config.mk:
#   -O0 -mp=gpu -stdpar=gpu -gpu=mem:separate,cc80,sm_80 -Mnofma -Minfo=all -r8
MRE_FLAGS = -O2 -mp=gpu -stdpar=gpu -gpu=mem:separate,cc80,sm_80 -Mnofma -r8

MRE_TYPES_OBJ = mre_ocean_types_mod.o

.PHONY: all clean mre1 mre1a mre2 mre_gradke mre_ptr_silent mre_ptr_silent_reshape mre_extptr mre_extptr_reshape

all: $(EXE)

# MRE 1: LLC "use of undefined value '%g'" compile-time error.
# Reproduces the bug in calc_slope_functions_using_just_e when G (a derived type
# with nested type + pointer members) is NOT mapped via !$omp target enter data
# but its allocatable members are accessed in an outer stdpar do concurrent.
# Expected: compilation FAILS with "error: use of undefined value '%g'" from llc.
mre1: $(MRE_TYPES_OBJ) mre1_llc_undefined_value.o
	$(FC) $(MRE_FLAGS) -o mre1.exe $(MRE_TYPES_OBJ) mre1_llc_undefined_value.o

# MRE 1a: Variant of MRE 1 without !$OMP parallel do, with outer do concurrent(j) loops restored.
# Matches calc_slope_functions_using_just_e after removing !$OMP parallel do and restoring the
# outer do concurrent(j=js:je) and do concurrent(J=js-1:je) loops (which were commented out and
# replaced with plain do j loops in the current repo). The k-loop remains a plain sequential loop.
# NOTE: With the simplified grid_type (8 members), this compiles successfully. The LLC error
# requires the full ocean_grid_type complexity (30+ allocatable members) to trigger without OMP.
# See the file header for a detailed explanation of why OMP is the essential ingredient.
mre1a: $(MRE_TYPES_OBJ) mre1a_llc_no_omp_parallel.o
	$(FC) $(MRE_FLAGS) -o mre1a.exe $(MRE_TYPES_OBJ) mre1a_llc_no_omp_parallel.o


# Reproduces the bug in calc_Visbeck_coeffs_old where local(H_u)/local(H_v)
# in an outer stdpar do concurrent causes a segfault when H_u/H_v are also
# explicitly mapped via !$omp target enter data.
# Expected: compilation succeeds; execution on a GPU node SEGFAULTS.
mre2: $(MRE_TYPES_OBJ) mre2_do_concurrent_local.o
	$(FC) $(MRE_FLAGS) -o mre2.exe $(MRE_TYPES_OBJ) mre2_do_concurrent_local.o

# MRE for gradKE failing to inline in MOM_CoriolisAdv.F90.
# Reproduces the 'array reshaping not enabled' inlining failure seen in MOM6.
#
# Without reshape (shows the failure, matching the MOM6 symptom):
#   make mre_gradke
# Expected: "subprogram not inlined -- array reshaping not enabled: gradke_mre2, argument 1"
#
# With reshape (demonstrates the fix):
#   make mre_gradke_reshape
# Expected: "gradke_mre2 inlined, size=209, ..."
mre_gradke: MRE_gradKE_failing_inline_MOM6.F90
	$(FC) $(MRE_FLAGS) -Minfo=inline -Minline=name:gradKE_mre2 \
	    MRE_gradKE_failing_inline_MOM6.F90 -o mre_gradke.exe

mre_gradke_reshape: MRE_gradKE_failing_inline_MOM6.F90
	$(FC) $(MRE_FLAGS) -Minfo=inline -Minline=name:gradKE_mre2,reshape \
	    MRE_gradKE_failing_inline_MOM6.F90 -o mre_gradke_reshape.exe

# MRE testing whether pointer members in a derived-type argument cause silent
# inlining failure in nvfortran.
#
# Tests 4 variants of gradKE:
#   A: G=no-ptr, CS=no-ptr  -> should inline  ("gradke_a inlined, size=...")
#   B: G=no-ptr, CS=with-ptr -> hypothesis: SILENT
#   C: G=with-ptr, CS=no-ptr -> hypothesis: SILENT
#   D: G=with-ptr, CS=with-ptr -> hypothesis: SILENT (mirrors real MOM6 gradKE)
#
# Without reshape (shows array reshaping failures for any variant that would inline):
#   make mre_ptr_silent
#
# With reshape (tests the pointer-member hypothesis):
#   make mre_ptr_silent_reshape
mre_ptr_silent: MRE_pointer_member_silent_noinline.F90
	$(FC) $(MRE_FLAGS) -Minfo=inline \
	    -Minline=name:gradKE_A,name:gradKE_B,name:gradKE_C,name:gradKE_D \
	    MRE_pointer_member_silent_noinline.F90 -o mre_ptr_silent.exe

mre_ptr_silent_reshape: MRE_pointer_member_silent_noinline.F90
	$(FC) $(MRE_FLAGS) -Minfo=inline \
	    -Minline=name:gradKE_A,name:gradKE_B,name:gradKE_C,name:gradKE_D,reshape \
	    MRE_pointer_member_silent_noinline.F90 -o mre_ptr_silent_reshape.exe

# ==========================================================================
# MRE: Silent non-inlining when derived type contains pointer to external type
#
# Root cause: when a subroutine argument's derived type contains a POINTER
# to any type defined in a SEPARATELY COMPILED MODULE, nvfortran silently
# skips inlining — no "inlined" or "subprogram not inlined" message appears.
#
# Two variants are compiled in a single executable:
#   Variant A (grid_no_ptr):   no external pointer → gradKE_no_ptr IS inlined
#   Variant B (grid_with_ptr): pointer to external_type → gradKE_with_ptr SILENT
#
# Expected -Minfo=inline output:
#   gradke_no_ptr:
#   gradke_with_ptr:
#   caller_no_ptr:
#     N, gradke_no_ptr inlined, size=..., ...    <- message present (inlined)
#   caller_with_ptr:
#                                                 <- NO message at all (silent)
#
# To build and see the difference:
#   make mre_extptr_reshape   (with -Minline=name:...,reshape; shows the contrast)
#   make mre_extptr           (without reshape; both fail but Variant A says why)
# ==========================================================================
external_type_mod.o: external_type_mod.F90
	$(FC) $(MRE_FLAGS) -c $<

mre_extptr_reshape: external_type_mod.o MRE_external_ptr_silent_noinline.F90
	$(FC) $(MRE_FLAGS) -Minfo=inline \
	    -Minline=name:gradKE_no_ptr,reshape \
	    -Minline=name:gradKE_with_ptr,reshape \
	    external_type_mod.o MRE_external_ptr_silent_noinline.F90 -o mre_extptr_reshape.exe

mre_extptr: external_type_mod.o MRE_external_ptr_silent_noinline.F90
	$(FC) $(MRE_FLAGS) -Minfo=inline \
	    -Minline=name:gradKE_no_ptr \
	    -Minline=name:gradKE_with_ptr \
	    external_type_mod.o MRE_external_ptr_silent_noinline.F90 -o mre_extptr.exe

$(MRE_TYPES_OBJ): mre_ocean_types_mod.F90
	$(FC) $(MRE_FLAGS) -c $<

mre1_llc_undefined_value.o: mre1_llc_undefined_value.F90 $(MRE_TYPES_OBJ)
	$(FC) $(MRE_FLAGS) -c $<

mre1a_llc_no_omp_parallel.o: mre1a_llc_no_omp_parallel.F90 $(MRE_TYPES_OBJ)
	$(FC) $(MRE_FLAGS) -c $<

mre2_do_concurrent_local.o: mre2_do_concurrent_local.F90 $(MRE_TYPES_OBJ)
	$(FC) $(MRE_FLAGS) -c $<

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJ) $(EXE) *.mod *.s *.optrpt
	rm -f $(MRE_TYPES_OBJ) mre1_llc_undefined_value.o mre1a_llc_no_omp_parallel.o mre2_do_concurrent_local.o
	rm -f mre1.exe mre1a.exe mre2.exe mre_gradke.exe mre_gradke_reshape.exe
	rm -f mre_ptr_silent.exe mre_ptr_silent_reshape.exe
	rm -f external_type_mod.o mre_extptr.exe mre_extptr_reshape.exe
