FC = nvfortran
FFLAGS = -O2 -mp=gpu -stdpar=gpu -gpu=mem:separate -Mnofma -Minfo=all -Minline

SRC = MRE_subroutine_call_in_teams_region.F90 MRE_subroutine_call_in_teams_region_program.F90
OBJ = $(SRC:.F90=.o)
EXE = subroutine_call_in_teams_region_from_makefile_O2.exe

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJ) $(EXE) *.mod *.s *.optrpt
