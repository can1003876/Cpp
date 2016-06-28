//============================================================================
// Name        : func_main.h
// Author      : 1003876
// Version     :
// Copyright   : Your copyright notice
// Description : General functions for TMatrix class
//============================================================================


#ifndef FUNC_MAIN_H_
#define FUNC_MAIN_H_

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

template <typename T>
TMatrix<T> operator*(T c, TMatrix<T>& rMatr){
	TMatrix<double> Amat(rMatr.getShape(0), rMatr.getShape(1), to_string(c) + "*" + rMatr.getName());
	for (int i=1;i<=Amat.getShape(0);i++){
		for (int j=1;j<=Amat.getShape(1);j++){
			Amat.setElement(i,j, c*rMatr.getElement(i,j) );
		}
	}
	return Amat;
}

template <typename T>
T max_value(const TMatrix<T>& rMatr){
	T temp = 0;
	for (int i=1; i<=rMatr.getShape(0); i++){
		for (int j=1; j<=rMatr.getShape(1); j++){
			if (temp < rMatr.getElement(i,j)){
				temp = rMatr.getElement(i,j);
			}
		}
	}
	return temp;
}

template <typename T>
double inner(const TMatrix<T>& vec1, const TMatrix<T>& vec2){
	assert(vec1.getShape(0) == vec2.getShape(0) && vec1.getShape(1) == 1 && vec2.getShape(1) == 1);
	double temp = 0;
	for (int i=1; i<=vec1.getShape(0); i++){
		temp += vec1.getElement(i,1)*vec2.getElement(i,1);
	}
	return temp;
}

template <typename T>
double norm(const TMatrix<T>& vec1, double p){
	assert(vec1.getShape(1) == 1);
	double temp = 0.0;
	for (int i=1; i<=vec1.getShape(0); i++){
		temp += (double) pow(vec1.getElement(i,1),p);
	}
	return pow(temp,1/p);
}

template <typename T>
double norm(const TMatrix<T>& vec1){
	assert(vec1.getShape(1) == 1);
	double temp = 0.0;
	for (int i=1; i<=vec1.getShape(0); i++){
		temp += (double) vec1.getElement(i,1)*vec1.getElement(i,1);
	}
	return sqrt(temp);
}

template <typename T>
T sum(const TMatrix<T>& rMatr){
	T temp = 0.0;
	for (int i=1; i<=rMatr.getShape(0); i++){
		for (int j=1; j<=rMatr.getShape(1); j++){
			temp += rMatr.getElement(i,j);
		}
	}
	return temp;
}

// either makes diagonal matrix or extracts diagonal from the matrix
template <typename T>
TMatrix<T> diag(const TMatrix<T>& rMatr){
	if (rMatr.getShape(0) == rMatr.getShape(1)){
		TMatrix<T> newMatrix(rMatr.getShape(0), 1,  "diag(" + rMatr.getName() + ")");
		for (int k=1; k<=newMatrix.getShape(0);k++){
			newMatrix(k,1) = rMatr.getElement(k,k);
		}
		return newMatrix;
	}else if (rMatr.getShape(0) == 1){
		TMatrix<T> newMatrix(rMatr.getShape(1), rMatr.getShape(1),"diag(" + rMatr.getName() + ")");

		for (int k=1; k<=newMatrix.getShape(0);k++){
			newMatrix(k,k) = rMatr.getElement(1,k);
		}
		return newMatrix;
	}else if (rMatr.getShape(1) == 1){
		TMatrix<T> newMatrix(rMatr.getShape(0), rMatr.getShape(0),"diag(" + rMatr.getName() + ")");

		for (int k=1; k<=newMatrix.getShape(1);k++){
			newMatrix(k,k) = rMatr.getElement(k,1);
		}
		return newMatrix;
	}else{
		cout << "Invalid input. Returning empty 1x1 matrix.";
		TMatrix<T> temp(1, 1, "empty 1x1 matrix");
		return temp;
	}
}

// This function returns new matrix that is transpose of the argument.
template <typename T>
TMatrix<T> transpose(const TMatrix<T>& rOther){
	TMatrix<T> rMatr(rOther.getShape(1), rOther.getShape(0),  rOther.getName() + "^T");
	for (int i=1; i<=rMatr.getShape(0); i++){
		for (int j=1; j<=rMatr.getShape(1); j++){
			rMatr.setElement(i,j,rOther.getElement(j,i));
		}
	}
	return rMatr;
}

// Check whether is symmetric.
template <typename T>
bool symmetric(const TMatrix<T>& rMatr){
	for (int i=1; i<=rMatr.getShape(0); i++){
		for (int j=1; j<=i-1; j++){
			if (rMatr.getElement(i,j) != rMatr.getElement(j,i)){
				return false;
			}
		}
	}
	return true;
}

// Returns Frobenius norm of how rMatr differs from a symmetric matrix.
template <typename T>
double symmetricF(const TMatrix<T>& rMatr){
	double temp=0;
	for (int i=1; i<=rMatr.getShape(0); i++){
		for (int j=1; j<=i-1; j++){
			temp += pow(2*abs(rMatr.getElement(i,j) - rMatr.getElement(j,i)),2);
		}
	}
	return sqrt(temp);
}

// Returns whether matrix is diagonally strictly dominant
template <typename T>
double diagSDom(const TMatrix<T>& rMatr){
	if (diagCDom(rMatr)){
		return true;
	}
	if (diagRDom(rMatr)){
		return true;
	}
	return false;
}

// Returns whether matrix is diagonally row dominant
template <typename T>
bool diagRDom(const TMatrix<T>& rMatr){
	double temp = 0;
	for (int i=1; i<=rMatr.getShape(0); i++){
		temp = 0;
		for (int j=1; j<=rMatr.getShape(1); j++){
			if (i != j){
				temp += abs(rMatr.getElement(i,j));
			}else{
				temp -= abs(rMatr.getElement(i,j));
			}
		}
		if (temp >= 0){
			return false;
		}
	}
	return true;
}

// Returns whether matrix is diagonally row dominant
template <typename T>
bool diagCDom(const TMatrix<T>& rMatr){
	double temp = 0;
	for (int j=1; j<=rMatr.getShape(1); j++){
		temp = 0;
		for (int i=1; i<=rMatr.getShape(0); i++){
			if (j != i){
				temp += abs(rMatr.getElement(i,j));
			}else{
				temp -= abs(rMatr.getElement(i,j));
			}
		}
		if (temp >= 0){
			return false;
		}
	}
	return true;
}

// Random matrix generator
// reference: http://www.cplusplus.com/reference/random/normal_distribution/normal_distirbution/
template <typename T>
TMatrix<T> RandomGauss(int sizeRows, int sizeColumns, double mean, double sd){
	TMatrix<T> rMatr(sizeRows, sizeColumns, "Random gauss matrix");
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	normal_distribution<double> distribution(mean, sd);
	for (int i=1; i<=rMatr.getShape(0); i++){
		for (int j=1; j<=rMatr.getShape(1); j++){
			rMatr(i,j) = (T) distribution(generator);
		}
	}
	return rMatr;
}

// Saves matrix into csv format. Load in MATLAB using csvread('A.csv');
template <typename T>
int saveMat(const TMatrix<T>& rMatr, string fileName){
	ofstream file;
	file.open(fileName);
	for (int i=1; i<=rMatr.getShape(0); i++){
		for (int j=1; j<=rMatr.getShape(1); j++){
			file << rMatr.getElement(i,j);
			if (j != rMatr.getShape(1)){
				file << ",";
			}
		}
		if (i != rMatr.getShape(0)){
			file << "\n";
		}
	}
	file.close();
	return 0;
}

// Saves matrix into csv format. File is named after the matrix.
// Load it in matlab using csvread('A.csv');
template <typename T>
int saveMat(const TMatrix<T>& rMatr){
	ofstream file;
	file.open(rMatr.getName()+".csv");
	for (int i=1; i<=rMatr.getShape(0); i++){
		for (int j=1; j<=rMatr.getShape(1); j++){
			file << rMatr.getElement(i,j);
			if (j != rMatr.getShape(1)){
				file << ",";
			}
		}
		if (i != rMatr.getShape(0)){
			file << "\n";
		}
	}
	file.close();
	return 0;
}

// We need to take copy of the matrix since we are going to alter it in LU.
template <typename T>
TMatrix<double> solveLU(TMatrix<T> Amat, const TMatrix<T>& bvec, bool pivoting){
	// Assert we have square matrix and vectors.
	assert(Amat.getShape(0) == bvec.getShape(0) && Amat.getShape(0) == Amat.getShape(1) && bvec.getShape(1) == 1);

	// Create copy of bvec, which will be changed to be the systems solution.
	TMatrix<double> xvec(bvec, "xvec");

	// Get Ux = Linv b
	for (int i=1; i<=Amat.getShape(1); i++){
		if (pivoting){
			cout << "Pivoting not implemented yet.";
		}
		for (int j=i+1; j<=Amat.getShape(0); j++){
			for (int k=i+1; k<=Amat.getShape(1); k++){
				Amat(j,k) -= Amat(j,i) / Amat(i,i) * Amat(i,k);
			}
			xvec(j,1) -= Amat(j,i) / Amat(i,i) * xvec(i,1);
			Amat(j,i) = 0;
		}
	}
	// Solve Ux = Linv b (Backsolve)
	for (int j=Amat.getShape(0); j>=1; j--){
		double temp = 0;
		for (int i=Amat.getShape(0); i>=j+1; i--){
			temp += Amat(j,i)*xvec(i,1);
		}
		xvec(j,1) = (xvec(j,1) - temp)/Amat(j,j);
	}
	return xvec;
}

// Reference: https://en.wikipedia.org/wiki/Gauss–Seidel_method#Algorithm
template <typename T>
TMatrix<double> solveGS(TMatrix<T> Amat, const TMatrix<T>& bvec, int max_iter){
	// Assert we have square matrix and vectors.
	assert(Amat.getShape(0) == bvec.getShape(0) && Amat.getShape(0) == Amat.getShape(1) && bvec.getShape(1) == 1);

	// Create copy of bvec, which will be changed to be the systems solution.
	TMatrix<double> xvec(bvec, "xvec");
	TMatrix<double> xvec_old(bvec.getShape(0), 1, "xvec_old");

	int iter = 1;
	do {
		for (int i=1; i<=Amat.getShape(0); i++){
			xvec_old = xvec;
			double temp = 0;
			for (int j=1; j<=Amat.getShape(1); j++){
				if (j != i){
					temp += Amat(i,j)*xvec(j,1);
				}
			}
			xvec(i,1) = (bvec.getElement(i,1)-temp)/Amat(i,i);
			iter++;
		}
	}while( norm(xvec-xvec_old,2.0) > 1.e-16 && iter < max_iter);
	return xvec;
}

template <typename T>
TMatrix<double> solveJAC(TMatrix<T> Amat, const TMatrix<T>& bvec, int max_iter){
	// Assert we have square matrix and vectors.
	assert(Amat.getShape(0) == bvec.getShape(0) && Amat.getShape(0) == Amat.getShape(1) && bvec.getShape(1) == 1);

	// Create copy of bvec, which will be changed to be the systems solution.
	TMatrix<double> xvec(bvec, "xvec");
	TMatrix<double> xvec_old(bvec.getShape(0), 1, "xvec_old");
	int iter = 1;
	do {
		for (int i=1; i<=Amat.getShape(0); i++){
			xvec_old = xvec;
			double temp = 0;
			for (int j=1; j<=Amat.getShape(1); j++){
				if (j != i){
					temp += Amat(i,j)*xvec_old(j,1);
				}
			}
			xvec(i,1) = (bvec.getElement(i,1)-temp)/Amat(i,i);
		}
		iter++;
	}while( (norm(xvec-xvec_old,2.0) > 1.e-16) && iter < max_iter );
	return xvec;
}

//Reference: https://en.wikipedia.org/wiki/Conjugate_gradient_method#Example_code_in_MATLAB_.2F_GNU_Octave
template <typename T>
TMatrix<T> solveCG(TMatrix<T> Amat, TMatrix<T> bvec, int max_iter = 20){
	// Assert we have square matrix and vectors.
	assert(Amat.getShape(0) == bvec.getShape(0) && Amat.getShape(0) == Amat.getShape(1) && bvec.getShape(1) == 1);

	// Create copy of bvec, which will be changed to be the systems solution.
	TMatrix<double> xvec(bvec, "xvec");
	TMatrix<double> xvec_old(bvec.getShape(0), 1, "xvec_old");

	TMatrix<double> rvec(bvec - Amat*xvec, "rvec");
	TMatrix<double> pvec(rvec, "pvec");

	TMatrix<double> Ap(bvec.getShape(0), 1, "Ap");

	double alpha;
	double rtr = inner(rvec,rvec);
	double rtrnew;

	int iter = 1;

	do {
		xvec_old = xvec;
		Ap = Amat*pvec;
		alpha = rtr/inner(pvec,Ap);
		xvec = xvec + alpha*pvec;
		rvec = rvec - alpha*Ap;
		rtrnew = inner(rvec,rvec);
		if (rtrnew <= 1.e-16){
			return xvec;
		}
		pvec = rvec + (rtrnew/rtr)*pvec;
		rtr = rtrnew;
		iter++;
	}while( norm(xvec-xvec_old,2.0) > 1.e-16 && iter < max_iter );
	return xvec;
}

template <typename T>
TMatrix<T> solveSD(TMatrix<T> Amat, TMatrix<T> bvec, int max_iter){
	// Assert we have square matrix and vectors.
	assert(Amat.getShape(0) == bvec.getShape(0) && Amat.getShape(0) == Amat.getShape(1) && bvec.getShape(1) == 1);

	// Create copy of bvec, which will be changed to be the systems solution.
	TMatrix<double> xvec(bvec, "xvec");
	TMatrix<double> xvec_old(bvec.getShape(0), 1, "xvec_old");

	TMatrix<double> rvec(bvec.getShape(0), 1, "rvec");

	double alpha;
	double rtr;

	int iter = 1;

	do {
		rvec = bvec - Amat*xvec;
		xvec_old = xvec;
		rtr = inner(rvec,rvec);
		alpha = rtr/inner(rvec,Amat*rvec);
		xvec = xvec + alpha*rvec;
		iter++;
	}while( norm(xvec-xvec_old,2.0) > 1.e-16 && iter < max_iter );
	return xvec;
}


// Solve linear system using QR decomposition.
// Here we use modified Gram-Schmidt for QR. Trefethen pg 73
// Reference: http://ocw.mit.edu/courses/mathematics/18-335j-introduction-to-numerical-methods-fall-2010/lecture-notes/MIT18_335JF10_lec10a_hand.pdf
template <typename T>
TMatrix<T> solveQR(TMatrix<T> Amat, TMatrix<T> bvec){
	// Assert we have square matrix and vectors.
	assert(Amat.getShape(0) == bvec.getShape(0) && Amat.getShape(0) == Amat.getShape(1) && bvec.getShape(1) == 1);

	// Create copy of bvec, which will be changed to be the systems solution.
	TMatrix<double> xvec(bvec, "xvec");
	TMatrix<double> Rmat(Amat.getShape(1), Amat.getShape(1), "Rmat");
	TMatrix<double> Qmat(Amat.getShape(0), Amat.getShape(1), "Qmat");

	// Get QRx = b
	for (int i=1; i<=Amat.getShape(1); i++) {
		Rmat.setElement(i,i, norm(Amat(i)));
		for (int j=1; j<=Amat.getShape(0); j++) {
			Qmat.setElement(j,i, Amat.getElement(j,i)/Rmat.getElement(i,i)); // column-wise loop
		}
		for (int j=i+1; j<=Amat.getShape(1); j++) {
			Rmat.setElement(i,j, inner(Qmat(i), Amat(i)));
			for (int k=1; k<=Amat.getShape(0);k++){
				Amat.setElement(k,j, Amat.getElement(k,j) - Rmat(i,j)*Qmat(k,j) );
			}
		}
	}
	// transform b
	bvec = transpose(Qmat)*bvec;

	// Solve Rx = Q' b (Backsolve)
	for (int j=Amat.getShape(0); j>=1; j--){
		double temp = 0;
		for (int i=Amat.getShape(0); i>=j+1; i--){
			temp += Amat(j,i)*xvec(i,1);
		}
		xvec(j,1) = (xvec(j,1) - temp)/Amat(j,j);
	}
	return xvec;
}

// General solve for linear equations which chooses the best method.
template <typename T>
TMatrix<T> linsolve(const TMatrix<T> Amat, const TMatrix<T>& bvec){
	assert(Amat.getShape(0) == bvec.getShape(0) && Amat.getShape(0) == Amat.getShape(1) && bvec.getShape(1) == 1);
	TMatrix<double> xvec(bvec, "xvec");
	if (diagSDom(Amat)){
		xvec = solveGS(Amat, bvec, Amat.getShape(0)*5);
	}else{
		xvec = solveLU(Amat, bvec, false);
	}
	return xvec;
}

// This function solves system Ax = b.
template <typename T>
TMatrix<T> operator/(const TMatrix<T>& bvec, const TMatrix<T> Amat){
	assert(Amat.getShape(0) == bvec.getShape(0) && Amat.getShape(0) == Amat.getShape(1) && bvec.getShape(1) == 1);
	TMatrix<double> xvec(bvec, "xvec");
	xvec = linsolve(Amat, bvec);
	return xvec;
}


#endif /* FUNC_MAIN_H_ */
