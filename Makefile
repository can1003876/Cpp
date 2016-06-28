all: MatLib

#For debugging
OPT=-g -Wall
#For optimistaion
#OPT=-O

#All objects (except main) come from cpp and hpp 
%.o:	%.cpp %.hpp
	g++ ${OPT} -c -o $@ $<
#use_vectors relies on objects which rely on headers
MatLib:	ProjectMatlib.cpp TMatrix.o
		g++ ${OPT} -o ProjectMatlib ProjectMatlib.cpp TMatrix.o


clean:
	rm -f *.o *~ ProjectMatlib