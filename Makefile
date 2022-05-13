PNP: PNP_main
	-@mv PNP_main PNP
	
include ${PHG_MAKEFILE_INC}

PNP_main: PNP_main.o

PNP_main.o: PNP_main.c PNP_build_solver.c PNP_coefficient.c PNP_more.c  \
            PNP_func.h PNP_quad.h PNP.h

clean: 
	-@rm -f *.o PNP *.vtk
