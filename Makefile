FC = nvfortran
FFLAGS = -O2 -mp=gpu -stdpar=gpu -gpu=mem:separate -Mnofma -Minfo=all -Minline=reshape -DUSE_DERIVED_TYPES
#FFLAGS = -O2 -mp=gpu -stdpar=gpu -gpu=mem:separate -Mnofma -Minfo=all -Minline=reshape

SRC = derived_types_mod.F90 MRE_subroutine_call_in_teams_region.F90 chain_mod4.F90 chain_mod3.F90 chain_mod2.F90 chain_mod1.F90 MRE_subroutine_call_in_teams_region_program.F90
OBJ = $(SRC:.F90=.o)
EXE = many_calls_with_derived_types_O2.exe

# MRE flags matching the MOM6 GPU build flags from ocean_only/config.mk:
#   -O0 -mp=gpu -stdpar=gpu -gpu=mem:separate,cc80,sm_80 -Mnofma -Minfo=all -r8
MRE_FLAGS = -O0 -mp=gpu -stdpar=gpu -gpu=mem:separate,cc80,sm_80 -Mnofma -Minfo=all -r8

MRE_TYPES_OBJ = mre_ocean_types_mod.o

.PHONY: all clean mre1 mre2

all: $(EXE)

# MRE 1: LLC "use of undefined value '%g'" compile-time error.
# Reproduces the bug in calc_slope_functions_using_just_e when G (a derived type
# with nested type + pointer members) is NOT mapped via !$omp target enter data
# but its allocatable members are accessed in an outer stdpar do concurrent.
# Expected: compilation FAILS with "error: use of undefined value '%g'" from llc.
mre1: $(MRE_TYPES_OBJ) mre1_llc_undefined_value.o
	$(FC) $(MRE_FLAGS) -o mre1.exe $(MRE_TYPES_OBJ) mre1_llc_undefined_value.o

# MRE 2: Runtime segfault with local(array) in outer do concurrent.
# Reproduces the bug in calc_Visbeck_coeffs_old where local(H_u)/local(H_v)
# in an outer stdpar do concurrent causes a segfault when H_u/H_v are also
# explicitly mapped via !$omp target enter data.
# Expected: compilation succeeds; execution on a GPU node SEGFAULTS.
mre2: $(MRE_TYPES_OBJ) mre2_do_concurrent_local.o
	$(FC) $(MRE_FLAGS) -o mre2.exe $(MRE_TYPES_OBJ) mre2_do_concurrent_local.o

$(MRE_TYPES_OBJ): mre_ocean_types_mod.F90
	$(FC) $(MRE_FLAGS) -c $<

mre1_llc_undefined_value.o: mre1_llc_undefined_value.F90 $(MRE_TYPES_OBJ)
	$(FC) $(MRE_FLAGS) -c $<

mre2_do_concurrent_local.o: mre2_do_concurrent_local.F90 $(MRE_TYPES_OBJ)
	$(FC) $(MRE_FLAGS) -c $<

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJ) $(EXE) *.mod *.s *.optrpt
	rm -f $(MRE_TYPES_OBJ) mre1_llc_undefined_value.o mre2_do_concurrent_local.o
	rm -f mre1.exe mre2.exe
