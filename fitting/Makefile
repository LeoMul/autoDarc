FC=gfortran
FFLAGS=-O3 -fbacktrace -fcheck=all
SRC= \
    amoeba.f \
    amotry.f90 \
    minimization.f90\
	kind_param.f90 \
    trapz.f90\
	factorials.f90 \
    dsygevlib.f \
    matrix.f90\
    norm.f90\
    main.f90 \


OBJ=${SRC:.f90=.o}
PROGRAM = main
PRG_OBJ = $(PROGRAM).o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

main: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm *.o *.mod 