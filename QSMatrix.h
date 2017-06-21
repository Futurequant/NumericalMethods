// Header file for Matrix Class 
// Quanttiative Finance 
// Updated: 23/2/17
// Matthew Morris 
#ifndef __QS_MATRIX_H
#define __QS_MATRIX_H
#include <vector>

template <typename T> class QSMatrix 
{
	private:
		std:vector<std::vector <T> > mat;

		// Represent the number of rows and colulms 
		// unisgned only allows positive values
		unsigned rows;
		unsigned columns;

	public:
		QSMatrix(unsigned _rows, unsigned _cols, const T&, _initial);
		QSMatrix(const QSMatrix<T>& rhs);
		virtual ~QSMatrix()

		// Operator overloading for "Standard" mathematical matrix ops
		QSMatrix<T>& operator=(const QSMatrix<T>& rhs);

		// Matrix mathematical operations
		QSMatrix<T> operator+(const QSMatrix<T>& rhs);
		QSMatrix<T>& operator+=(const QSMatrix<T>& rhs);
		QSMatrix<T> operator-(const QSMatrix<T>& rhs);
		QSMatrix<T>& operator-=(const QSMatrix<T> rhs);
		QSMatrix<T> operator*(const QSMatrix<T>& rhs);
		QSMatrix<T>& operator*=(const QSMatrix<t> rhs);
		QSMatrix<T> transpose();

		// Matrix/Scalar Operations 
		QSMatrix<T> operator+(const T& rhs);
		QSMatrix<T> operator-(const T& rhs);
		QSMatrix<T> operator*(const T& rhs);
		QSMatrix<T> operator/(const T& rhs);

		// Matrix/Vector Operations 
		std::vector<T> operator*(const std::vector<T>& rhs);
		std::vector<T> diag_vec();

		// Access individual elements
		unsigned get_rows() const;
		unsigned get_cols() const;
};

// Include the matrix implementatin file for c++ template type deduction 
#include "matrix.cpp"

#endif