/*
 * TMatrix.h
 *
 *  Created on: 11 May 2016
 *      Author: 1003876
 */

#ifndef TMATRIX_H_
#define TMATRIX_H_

#include <cassert>
#include <string>
#include <iostream>
#include <random>
#include <math.h>

using namespace std;

template <typename T>
class TMatrix;

template<typename T>
ostream& operator<<(ostream& output, const TMatrix<T>& rMatr);

template <typename T>
class TMatrix {
private:
	T *mData;	// Matrix data (instead of: T **mData)
	int mShape[2];	// Dimensions of the matrix
	int mSize;	// Number of elements in the matrix
	string mName;	// Matrix name
public:
	// Functions:
	TMatrix(int sizeRows, int sizeColumns, string Name);	// Constructor
	TMatrix(int sizeRows, int sizeColumns, T elements, string Name);	// Constructor with specified elements;
	TMatrix(const TMatrix<T>& rOther);						// Copy Constructor
	TMatrix(const TMatrix<T>& rOther, string Name);			// Copy Constructor with specified name
	virtual ~TMatrix();		// Destructor

	// Getting/Changing data
	T& operator()(int i, int j);
	TMatrix<T> operator()(int ind);					// if only one index is specified, returns column/row vector.
	void swap(int i1, int j1, int i2, int j2);		// Swap two values in the matrix
	void setElement(int i, int j, T value);
	const T getElement(int i, int j) const;
	void transpose();

	// Overloaded explicit casts for type conversion
	TMatrix<int> toInt();
	TMatrix<double> toDouble();

	// Matrix operations
	TMatrix<T>& operator=(const TMatrix<T>& rMatr);
	bool operator==(const TMatrix<T>& rMatr);
	TMatrix<T> operator+(const TMatrix<T>& rOther);
	TMatrix<T> operator-(const TMatrix<T>& rOther);
	TMatrix<T> operator-();
	TMatrix<T> operator*(const TMatrix<T>& rOther);
	TMatrix<T> operator*(T scalar);
	//	TMatrix<T> operator/(const TMatrix<T>& rOther);
	void rMAdd(const TMatrix<T>& rOther);

	// Getting/Changing information
	void setName(string newName);
	const string getName() const;
	const int getShape(int dim) const;
	const int getShape() const;
	const int length() const;

	friend ostream& operator << <> (ostream& output, const TMatrix& rMatr);
};


// TMatrix implementation

// TMatrix constructor - initialize full of zeroes.
template <typename T>
TMatrix<T>::TMatrix(int sizeRows, int sizeColumns, string Name) {
	// Define size of a matrix
	mSize = sizeRows*sizeColumns;
	mShape[0] = sizeRows;
	mShape[1] = sizeColumns;
	// Matrix name - mostly used when printing
	mName = Name;

	// Matrix data - initialize and fill with zeroes.
	mData = new T[mSize];
	for (int i=0; i<mSize; i++)
	{
		mData[i] = (T) 0.0;
	}
}

// TMatrix constructor - initialize full of %elements.
template <typename T>
TMatrix<T>::TMatrix(int sizeRows, int sizeColumns, T elements, string Name) {
	// Define size of a matrix
	mSize = sizeRows*sizeColumns;
	mShape[0] = sizeRows;
	mShape[1] = sizeColumns;
	// Matrix name - mostly used when printing
	mName = Name;

	// Matrix data - initialize and fill with zeroes.
	mData = new T[mSize];
	for (int i=0; i<mSize; i++)
	{
		mData[i] = (T) elements;
	}
}

// Copy Constructor of TMatrix
template <typename T>
TMatrix<T>::TMatrix(const TMatrix<T>& rOther) {
	  mSize = rOther.mSize;
	  mShape[0] = rOther.mShape[0];
	  mShape[1] = rOther.mShape[1];
	  mName = rOther.mName + "-copy";
	  // cout << "Name is " << mName << "\n"; // This is for checking the names
	  mData = new T[mSize];
	  for (int i=0; i<rOther.mSize; i++){
		  mData[i] = rOther.mData[i];
	  }
}

// Copy Constructor of TMatrix with specified name
template <typename T>
TMatrix<T>::TMatrix(const TMatrix<T>& rOther, string Name) {
	  mSize = rOther.mSize;
	  mShape[0] = rOther.mShape[0];
	  mShape[1] = rOther.mShape[1];
	  mName = Name;
	  // cout << "Name is " << mName << "\n"; // This is for checking the names
	  mData = new T[mSize];
	  for (int i=0; i<rOther.mSize; i++){
		  mData[i] = rOther.mData[i];
	  }
}

// TMatrix destructor
template <typename T>
TMatrix<T>::~TMatrix() {
	delete[] mData;
}

template <typename T>
void TMatrix<T>::setElement(int i, int j, T value){
	assert(j>0 && i>0);
	mData[(j-1) + (i-1)*mShape[1]] = value;
}

// This function returns reference(constant) to mData. Can be only used for reading mData.
template <typename T>
const T TMatrix<T>::getElement(int i, int j) const{
	assert(j>0 && i>0);
	return mData[(j-1) + (i-1)*mShape[1]];
}

// Swaps two values in a matrix
template <typename T>
void TMatrix<T>::swap(int i1, int j1, int i2, int j2){
	mData[(j1-1) + (i1-1)*mShape[1]] += mData[(j2-1) + (i2-1)*mShape[1]];
	mData[(j2-1) + (i2-1)*mShape[1]]  = mData[(j1-1) + (i1-1)*mShape[1]] - mData[(j2-1) + (i2-1)*mShape[1]];
	mData[(j1-1) + (i1-1)*mShape[1]] -= mData[(j2-1) + (i2-1)*mShape[1]];
}

// This function returns reference(constant) to mData. Can be only used for reading mData.
template <typename T>
void TMatrix<T>::transpose(){
	for (int i=1; i<=getShape(0); i++){
		for (int j=1; j<=i-1; j++){
			swap(i,j,j,i);
		}
	}
}

// This function returns reference to mData. We can read or set values using this.
template <typename T>
T& TMatrix<T>::operator()(int i, int j){
	assert(j>0 && i>0);
	return mData[(j-1) + (i-1)*mShape[1]];
}

// This function returns new TMatrix object which corresponds either to the ind-column/row (based on sign).
template <typename T>
TMatrix<T> TMatrix<T>::operator()(int ind){
	if (ind > 0){
		// pick ind-th column vector
		string newName = mName + "-" + to_string(ind)+"column";
		TMatrix<T> temp(mShape[0], 1, newName);
		for (int j=1; j<=mShape[0]; j++){
			temp(j,1) = (*this)(j,ind);
		}
		return temp;
	}else if (ind < 0){
		// pick ind-th row
		string newName = mName + "-" + to_string(-ind)+"row";
		TMatrix<T> temp(1, mShape[1], newName);
		for (int j=1; j<=mShape[1]; j++){
			temp(1,j) = (*this)(-ind,j);
		}
		return temp;
	}else{
		cout << "Invalid index. Returning empty 1x1 matrix.";
		TMatrix<T> temp(1, 1, "empty 1x1 matrix");
		return temp;
	}
}

// This functions returns copy of the matrix with int mData type.
template <typename T>
TMatrix<int> TMatrix<T>::toInt() {
	TMatrix<int> rMatr(mShape[0], mShape[1], mName + "-int");
	for (int i=1; i<=mShape[0]; i++){
			for (int j=1; j<=mShape[1]; j++){
				rMatr.setElement(i,j, round( getElement(i,j) ) );
			}
	}
	return rMatr;
}

// This functions returns copy of the matrix with double mData type.
template <typename T>
TMatrix<double> TMatrix<T>::toDouble() {
	TMatrix<double> rMatr(mShape[0], mShape[1], mName + "-double");
	for (int i=1; i<=mShape[0]; i++){
			for (int j=1; j<=mShape[1]; j++){
				rMatr.setElement(i,j, (double) getElement(i,j));
			}
	}
	return rMatr;
}

// This function creates a copy of the matrix rMatr.
template <typename T>
TMatrix<T>& TMatrix<T>::operator=(const TMatrix<T>& rMatr){
	for (int i=0; i<rMatr.mSize; i++){
		mData[i] = rMatr.mData[i];
	}
	return *this;
}

// Returns true if all elements are same
template <typename T>
bool TMatrix<T>::operator==(const TMatrix<T>& rMatr){
	bool same = true;
	for (int i=0; i<rMatr.mSize; i++){
		if (mData[i] != rMatr.mData[i]){
			same = false;
		};
	}
	return same;
}

// Returns new TMatrix object that is sum of two matrices.
template <typename T>
TMatrix<T> TMatrix<T>::operator+(const TMatrix<T>& rOther){
	assert(mShape[1] == rOther.mShape[1] && mShape[0] == rOther.mShape[0]);
	TMatrix<T> rMatr(mShape[0], mShape[1],  mName + "+" + rOther.getName());
	for (int i=0; i< mSize; i++){
		rMatr.mData[i] = mData[i] + rOther.mData[i];
	}
	return rMatr;
}

// This function returns new TMatrix object that is subtraction of two matrices.
template <typename T>
TMatrix<T> TMatrix<T>::operator-(const TMatrix<T>& rOther){
	TMatrix<T> rMatr(mShape[0], mShape[1],  mName + "-" + rOther.getName());
	for (int i=0; i< mSize; i++){
		rMatr.mData[i] = mData[i] - rOther.mData[i];
	}
	return rMatr;
}

// This function returns new TMatrix object that is negative matrix.
template <typename T>
TMatrix<T> TMatrix<T>::operator-(){
	TMatrix<T> rMatr(mShape[0], mShape[1], "-" + mName);
	for (int i=0; i< mSize; i++){
		rMatr.mData[i] = - mData[i];
	}
	return rMatr;
}

// This function returns new TMatrix object that is matrix multiplication of two matrices.
template <typename T>
TMatrix<T> TMatrix<T>::operator*(const TMatrix<T>& rOther){
	assert(mShape[1] == rOther.mShape[0]);
	TMatrix rMatr(mShape[0], rOther.mShape[1], mName + "*" + rOther.getName());
	for (int i=1; i<=rMatr.mShape[0]; i++){
		for (int j=1; j<=rMatr.mShape[1]; j++){
			T temp=0;
			for (int k=1; k<=mShape[1]; k++){
				temp = temp + getElement(i,k)*rOther.getElement(k,j);
			}
			rMatr.setElement(i,j,temp);
		}
	}
	return rMatr;
}

// This function returns new TMatrix object that is rMatr=scalar*matrix.
template <typename T>
TMatrix<T> TMatrix<T>::operator*(T scalar){
	TMatrix rMatr(mShape[0], mShape[1], to_string(scalar) + "*" + mName);
	for (int i=0; i< mSize; i++){
		rMatr.mData[i] =  scalar*mData[i];
	}
	return rMatr;
}

// This is right-multiplication.
// Takes reference of other matrix and multiplies with it from the right.
// Original matrix is changed, but no copying is happening, so it is faster than operator*.
template <typename T>
void TMatrix<T>::rMAdd(const TMatrix<T>& rOther){
	assert(mShape[1] == rOther.mShape[1] && mShape[0] == rOther.mShape[0]);
	for (int i=0; i< mSize; i++){
		mData[i] = mData[i] + rOther.mData[i];
	}
}

template <typename T>
const string TMatrix<T>::getName() const{
	return mName;
}

template <typename T>
void TMatrix<T>::setName(string newName){
	mName = newName;
}

template <typename T>
const int TMatrix<T>::getShape(int dim) const{
	return mShape[dim];
}

template <typename T>
const int TMatrix<T>::getShape() const{
	return mSize;
}

template <typename T>
const int TMatrix<T>::length() const{
	return mSize;
}

// Takes reference to a matrix and prints its data.
template <typename T>
ostream& operator<<(ostream& output, const TMatrix<T>& rMatr){
	output << rMatr.getName() << " = \n";
	for (int i=1; i<=rMatr.mShape[0]; i++){
		output << "(" << rMatr.getElement(i,1);
		for (int j=2; j<=rMatr.mShape[1]; j++){
			output << ", " << rMatr.getElement(i,j);
		}
		output << ") \n";
	}
	return output;
}



#endif /* TMATRIX_H_ */
