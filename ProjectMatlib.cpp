//============================================================
// Name        :ProjectMatlib.cpp
// Author      :1003876
// Version     :
// Copyright   :Your copyright notice
// Description :Showcase of TMatrix class + solver comparisons
//============================================================

#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <string>
#include <iostream>
#include <math.h>
#include <time.h>

#include "TMatrix.h"
#include "func_main.h"


using namespace std;

int main(void) {


// Code listing 1 (lines 26 to 51):
// 2x3 double matrix "A":
TMatrix<double> A(RandomGauss<double>(2, 3, 1.0, 3.0), "A");
A(2,2) = 1; // Change A_{2,2} == 1:

// Create integer matrix B by rounding elements of A:
TMatrix<int> B(transpose(A).toInt(), "B");
cout << A;
//A =
//(0.115494, 4.04615, 7.23705)
//(-0.752548, 1, 0.569711)

cout << B;
//B =
//(0, -1)
//(4, 1)
//(7, 1)

// Transpose B, cast it as double and add it to A:
cout << A+transpose(B).toDouble();
//A+B^T-double =
//(0.115494, 8.04615, 14.237)
//(-1.75255, 2, 1.56971)

cout << A*B.toDouble();		//name : "A*B-double"
cout << A.toInt()*B;		//name : "A-int*B"

// Code listing 2 (lines 53 to 62):
cout << A;
//A =
//(0.115494, 4.04615, 7.23705)
//(-0.752548, 1, 0.569711)

cout << A(-1);
//A-1row = (0.115494, 4.04615, 7.23705)
cout << transpose(B(1).toInt());	// B = A^T.toInt()
//B-1column^T = (0, 4, 7)

// Code listing 3 (lines 64 to 72):
TMatrix<double> C(RandomGauss<double>(5e3, 1e3, 0.0, 3.0), "C");
TMatrix<double> D(RandomGauss<double>(5e3, 1e3, 0.0, 3.0), "D");
clock_t start = clock ();
C = C+D;	// It took: 0.317092 sec.
printf("It took: %Lf sec.\n", (long double) (clock() - start)/ CLOCKS_PER_SEC );
start = clock ();
C.rMAdd(D);	//It took: 0.085300 sec.
printf("It took: %Lf sec.\n", (long double) (clock() - start)/ CLOCKS_PER_SEC );

// Code listing 4 (lines 74 to 89):
TMatrix<int> E(2, 2, 1, "E");
cout << E;
//E =
//(1, 1)
//(1, 1)
cout << transpose(diag(E));
//diag(E)^T = (1, 1)
cout << diag( diag(E) );
//diag(diag(E)) =
//(1, 0)
//(0, 1)
saveMat(diag( diag(E) + transpose(B(-1)) ), "diag.csv");
// diag.csv can be loaded in MATLAB with: csvread('diag.csv');
cout << sum(E) + norm(diag(E)) - sqrt(2); cout << "\n";
// 4

// Code listing 6 (lines 91 to 98):
TMatrix<double> Amat(RandomGauss<double>(5, 5, 2.0, 3.0), "Amat");
TMatrix<double> xsol(RandomGauss<int>(5, 1, 3.0, 4.0).toDouble(), "xsol");
TMatrix<double> bvec(Amat*xsol, "bvec");
cout << transpose(bvec/Amat);	// calls linsolve(Amat, bvec)
//xvec^T = (10, 4, -4, 5, 2)
cout << transpose(xsol);
//xsol^T = (10, 4, -4, 5, 2)

// Code listing 7 (lines 100 to 168):
// Start/end size, and step in iterations
int aN = 50, bN = 1000, dN = 50;
int MAX_ITER = 5000;
// Matrices to save times and errors
TMatrix<double> times(6,(int) (bN-aN)/dN, "times");
TMatrix<double> errors(6,(int) (bN-aN)/dN, "errors");

int k = 1; long double time = 0;
for (int N=aN; N<=bN; N=N+dN){
		TMatrix<double> T(RandomGauss<double>(N, N, 0.0, 1.0), "T");
		TMatrix<double> S(N, N, 3*N, "S");
		S = diag(diag(S));	// diagonal matrix full of 3*N elements
		T = T + transpose(T) + S;	// positive definite-diagonally dominant matrix
		TMatrix<double> x(RandomGauss<double>(N, 1, 2.0, 5.0).toDouble(), "x");
		TMatrix<double> b(T*x, "b");

//		LU solver
		start = clock ();
		TMatrix<double> xLU(solveLU(T, b, false), "xLU");
		time = (long double) (clock() - start)/ CLOCKS_PER_SEC;
		printf("Finished solveLU. L2 error is: %1.3e and it took: %Lf sec. \n", norm(x-xLU), time);
		times(1,k) = time;
		errors(1,k) = norm(x-xLU);

//		QR solver
		start = clock ();
		TMatrix<double> xQR(solveQR(T, b), "xQR");
		time = (long double) (clock() - start)/ CLOCKS_PER_SEC;
		printf("Finished solveQR. L2 error is: %1.3e and it took: %Lf sec. \n", norm(x-xQR), time);
		times(2,k) = time;
		errors(2,k) = norm(x-xQR);

//		JACobi solver
		start = clock ();
		TMatrix<double> xJAC(solveJAC(T, b, MAX_ITER), "xJAC");
		time = (long double) (clock() - start)/ CLOCKS_PER_SEC;
		printf("Finished solveJAC. L2 error is: %1.3e and it took: %Lf sec. \n", norm(x-xJAC), time);
		times(3,k) = time;
		errors(3,k) = norm(x-xJAC);

//		GaussSeidel solver
		start = clock ();
		TMatrix<double> xGS(solveGS(T, b, MAX_ITER), "xGS");
		time = (long double) (clock() - start)/ CLOCKS_PER_SEC;
		printf("Finished solveGS. L2 error is: %1.3e and it took: %Lf sec. \n", norm(x-xGS), time);
		times(4,k) = time;
		errors(4,k) = norm(x-xGS);

//		Steepest descent solver
		start = clock ();
		TMatrix<double> xSD(solveSD(T, b, MAX_ITER), "xSD");
		time = (long double) (clock() - start)/ CLOCKS_PER_SEC;
		printf("Finished solveSD. L2 error is: %1.3e and it took: %Lf sec. \n", norm(x-xSD), time);
		times(5,k) = time;
		errors(5,k) = norm(x-xSD);

//		Conjugate gradient solver
		start = clock ();
		TMatrix<double> xCG(solveCG(T, b, MAX_ITER), "xCG");
		time = (long double) (clock() - start)/ CLOCKS_PER_SEC;
		printf("Finished solveCG. L2 error is: %1.3e and it took: %Lf sec. \n", norm(x-xCG), time);
		times(6,k) = time;
		errors(6,k) = norm(x-xCG);

		k++;
	}
saveMat(times);
saveMat(errors);

printf("Fin !\n");
return 0;
}
