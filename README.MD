settings used to run MRE: mpifort -g -O0 -mp=gpu -stdpar=gpu -gpu=mem:separate -Mnofma -Minfo=all -r8 MRE_subroutine_call_in_teams_region.F90 -o subroutine_mre.exe
