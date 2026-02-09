FC = nvfortran
FFLAGS = -O2 -mp=gpu -stdpar=gpu -gpu=mem:separate -Mnofma -Minfo=all -Minline=reshape -DUSE_DERIVED_TYPES
#FFLAGS = -O2 -mp=gpu -stdpar=gpu -gpu=mem:separate -Mnofma -Minfo=all -Minline=reshape

SRC = derived_types_mod.F90 MRE_subroutine_call_in_teams_region.F90 chain_mod4.F90 chain_mod3.F90 chain_mod2.F90 chain_mod1.F90 MRE_subroutine_call_in_teams_region_program.F90
OBJ = $(SRC:.F90=.o)
EXE = many_calls_with_derived_types_O2.exe

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJ) $(EXE) *.mod *.s *.optrpt
