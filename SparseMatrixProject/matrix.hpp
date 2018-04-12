#pragma once
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <cmath>
#include <fstream>
#include <sstream>
#include <sstream>

/*
template <class T>
class matrix2D
{
	private:
		std::vector<T> data;
		int columnCount;
		int rowCount;

	public:
		void print() const;
		// int getColumnCount(return columnCount;)
		void print(std::ostream & os) const;
		T &operator()(int x, int y);
		friend std::ostream& operator<< (std::ostream& stream, const matrix2D<T>& matrix);
		friend matrix2D transpose(matrix2D<T>& matrix);
		friend matrix2D multiply(matrix2D<T> & m1, matrix2D<T> & m2);
		T* getData(int rowIdx);
		std::vector<T> getVector();
		matrix2D(int x, int y);
};
template <class T> class matrix2D;
template <class T>
std::ostream& operator<<(std::ostream& stream, const matrix2D<T>& matrix);
template <class T>
matrix2D<T> transpose(matrix2D<T>& matrix);
template <class T>
matrix2D<T> multiply(matrix2D<T> & m1, matrix2D<T> & m2);
*/

template <class T> class matrix2D;

template <class T>
std::ostream& operator<<(std::ostream& stream, const matrix2D<T>& matrix);

template <class T>
matrix2D<T> transpose(matrix2D<T>& matrix);

template <class T>
bool isValidToMultiply(matrix2D<T> & m1, matrix2D<T> & m2);

template <class T>
bool isValidToAdd(matrix2D<T> & m1, matrix2D<T> & m2);

template <class T>
matrix2D<T> multiply(matrix2D<T> & m1, matrix2D<T> & m2);

template <class T>
matrix2D<T> add(matrix2D<T> & m1, matrix2D<T> & m2);

template <class T>
class matrix2D
{
public:
	matrix2D();
	matrix2D(int rowCount, int colCount);
	matrix2D(int rowCount, int colCount, T fill);
	matrix2D(const std::string & mtxFile);
	void print() const;
	void print(std::ostream & os) const;
	T &operator()(int x, int y);
	T* getData(int rowIdx);
	std::vector<T> getVector();
	int getRowCount();
	int getColumnCount();
	int getNumElements();
	matrix2D<T> createSubMatrix(int startRow, int endRow, int startCol, int endCol);
	matrix2D<T> & getSubMatrix(int startRow, int endRow, int startCol, int endCol);
	void copyValues(matrix2D<T> & from, int startRow, int endRow);
	void copyValues(matrix2D<T> & from, int startRow, int endRow, int startCol, int endCol);
	void fillRandomly(const int min = 0, const int max = 5);
	void fillRandomlyLowerTriangular(const int min = 0, const int max = 5);
	void makePositiveDefinite();
	friend std::ostream& operator<< <>(std::ostream& stream, const matrix2D<T>& matrix);
	friend matrix2D<T> transpose <>(matrix2D<T>& matrix);
	friend bool isValidToMultiply <>(matrix2D<T> & m1, matrix2D<T> & m2);
	friend bool isValidToAdd <>(matrix2D<T> & m1, matrix2D<T> & m2);
	friend matrix2D<T> multiply <>(matrix2D<T> & m1, matrix2D<T> & m2);
	friend matrix2D<T> add <>(matrix2D<T> & m1, matrix2D<T> & m2);
private:
	std::vector<T> data;
	int columnCount;
	int rowCount;
	void resize(const unsigned int newRowCount, const unsigned int newColumnCount);
};
template <class T>
matrix2D<T>::matrix2D() : data(), rowCount(0), columnCount(0)
{

}

template <class T>
matrix2D<T>::matrix2D(int x, int y) : data(x*y), rowCount(x), columnCount(y)
{

}

template <class T>
matrix2D<T>::matrix2D(int x, int y, T fill) : data(x*y, fill), rowCount(x), columnCount(y)
{

}

template <class T>
matrix2D<T>::matrix2D(const std::string & mtxFile)
{
	std::ifstream file(mtxFile);
	std::string line;
	bool isHeaderLine = true;
	while (std::getline(file, line))
	{
		if (line.empty() || line[0] == '%')
		{
			continue;
		}
		std::istringstream iss(line);
		if (isHeaderLine)
		{
			int nNonZeros;
			iss >> this->rowCount;
			iss >> this->columnCount;
			iss >> nNonZeros;
			this->resize(rowCount, columnCount);
			isHeaderLine = false;
		}
		else
		{
			int row, col;
			T val;
			iss >> row;
			iss >> col;
			iss >> val;
			row--;
			col--;
			(*this)(row, col) = val;
		}
	}
}

template <class T>
void matrix2D<T>::print() const
{
	for (int i = 0; i < rowCount; i++)
	{
		for (int j = 0; j < columnCount; j++)
		{
			std::cout << data[i * columnCount + j] << ", ";
		}
		std::cout << std::endl;
	}
}

template <class T>
void matrix2D<T>::print(std::ostream & os) const
{
	for (int i = 0; i < rowCount; i++)
	{
		for (int j = 0; j < columnCount; j++)
		{
			os << data[i * columnCount + j] << ", ";
		}
		os << std::endl;
	}
}

template <class T>
T& matrix2D<T>::operator()(int x, int y)
{
	return data[x * columnCount + y];
}

template <class T>
T* matrix2D<T>::getData(int rowIdx)
{
	return &data[rowIdx * columnCount];
}

template <class T>
std::vector<T> matrix2D<T>::getVector()
{
	return data;
}

template <class T>
int matrix2D<T>::getRowCount()
{
	return rowCount;
}

template <class T>
int matrix2D<T>::getColumnCount()
{
	return columnCount;
}

template <class T>
int matrix2D<T>::getNumElements()
{
	return rowCount * columnCount;
}

template <class T>
matrix2D<T> matrix2D<T>::createSubMatrix(int startRow, int endRow, int startCol, int endCol)
{
	if (startRow < 0 || startRow > rowCount - 1 || startCol < 0 || startCol > columnCount - 1
		|| endRow < 0 || endRow > rowCount - 1 || endCol < 0 || endCol > columnCount - 1
		|| startRow > endRow || startCol > endCol)
	{
		std::cerr << "createSubMatrix(). Invalid input. startRow=" << startRow << " endRow=" << endRow
			<< " startCol=" << startCol << " endCol=" << endCol << " this->rowCount=" << rowCount << " this->columnCount=" << columnCount << std::endl;
		return matrix2D<T>();
	}
	matrix2D<T> sub((endRow - startRow) + 1, (endCol - startCol) + 1);
	int row = 0, col = 0;
	for (int i = startRow; i <= endRow; i++)
	{
		for (int j = startCol; j <= endCol; j++)
		{
			sub(row, col) = (*(this))(i, j);
			col++;
		}
		col = 0;
		row++;
	}
	return sub;
}

template <class T>
matrix2D<T> & matrix2D<T>::getSubMatrix(int startRow, int endRow, int startCol, int endCol)
{
	if (startRow < 0 || startRow > rowCount - 1 || startCol < 0 || startCol > columnCount - 1
		|| endRow < 0 || endRow > rowCount - 1 || endCol < 0 || endCol > columnCount - 1
		|| startRow > endRow || startCol > endCol)
	{
		std::cerr << "createSubMatrix(). Invalid input. startRow=" << startRow << " endRow=" << endRow
			<< " startCol=" << startCol << " endCol=" << endCol << " this->rowCount=" << rowCount << " this->columnCount=" << columnCount << std::endl;
		return *this;
	}
	int newRowCount = (endRow - startRow) + 1;
	int newColumnCount = (endCol - startCol) + 1;
	int row = 0, col = 0;
	for (int i = startRow; i <= endRow; i++)
	{
		for (int j = startCol; j <= endCol; j++)
		{
			(*(this))(row, col) = (*(this))(i, j);
			col++;
		}
		col = 0;
		row++;
	}

	this->resize(newRowCount, newColumnCount);
	return *this;
}

template <class T>
void  matrix2D<T>::copyValues(matrix2D<T> & from, int startRow, int endRow)
{
	int fromRowCount = from.getRowCount();
	int fromColCount = from.getColumnCount();
	int numRowsToFillIn = (endRow - startRow) + 1;
	if (startRow < 0 || startRow > rowCount - 1 || endRow < 0 || endRow > rowCount - 1 ||
		fromColCount > columnCount)
	{
		std::cerr << "copyValues(). Invalid input. startRow=" << startRow << " endRow=" << endRow
			<< " this->rowCount=" << rowCount << " this->columnCount=" << columnCount << std::endl;
		return;
	}
	int row = 0, col = 0;
	for (int i = startRow; i <= endRow; i++)
	{
		for (int j = 0; j < fromColCount; j++)
		{
			(*(this))(i, j) = from(row, col);
			col++;
		}
		col = 0;
		row++;
		// Stop if we want to copy more rows than from matrix actually has
		if (row >= fromRowCount)
		{
			break;
		}
	}
}

template <class T>
void matrix2D<T>::copyValues(matrix2D<T> & from, int startRow, int endRow, int startCol, int endCol)
{
	// int numRowsToFillIn = (endRow - startRow) + 1;
	// int numColsToFillIn = (endCol - startCol) + 1;
	int fromRowCount = from.getRowCount();
	int fromColCount = from.getColumnCount();
	// if from is larger than "this" matrix then will fill in as many values as possible from top right
	if (startRow < 0 || startRow > rowCount - 1 || startCol < 0 || startCol > columnCount - 1
		|| endRow < 0 || endRow > rowCount - 1 || endCol < 0 || endCol > columnCount - 1
		|| startRow > endRow || startCol > endCol /* || from.getColumnCount() > columnCount
		|| from.getRowCount() > rowCount */)
	{
		std::cerr << "copyValues(). Invalid input. startRow=" << startRow << " endRow=" << endRow
			<< " startCol=" << startCol << " endCol=" << endCol << " this->rowCount=" << rowCount << " this->columnCount=" << columnCount << std::endl;
		return;
	}
	int row = 0, col = 0;
	for (int i = startRow; i <= endRow; i++)
	{
		for (int j = startCol; j <= endCol; j++)
		{
			(*(this))(i, j) = from(row, col);
			col++;
			// Stop if we want to copy more cols than from matrix actually has
			if (col >= fromColCount)
			{
				break;
			}
		}
		col = 0;
		row++;
		// Stop if we want to copy more rows than from matrix actually has
		if (row >= fromRowCount)
		{
			break;
		}
	}
}

template <class T>
void matrix2D<T>::fillRandomly(const int min, const int max)
{
	for (int i = 0; i < rowCount; i++)
	{
		for (int j = 0; j < columnCount; j++)
		{
			T randNum = (T)(min + (rand() % static_cast<int>(max - min + 1)));
			randNum = (randNum == 0.0) ? 1.0 : randNum;
			(*this)(i, j) = randNum;
		}
	}
}

template <class T>
void matrix2D<T>::fillRandomlyLowerTriangular(const int min, const int max)
{
	for (int i = 0; i < rowCount; i++)
	{
		for (int j = 0; j < columnCount; j++)
		{
			if (i >= j)
			{
				T randNum = (T)(min + (rand() % static_cast<int>(max - min + 1)));
				randNum = (randNum == 0.0) ? 1.0 : randNum;
				(*this)(i, j) = randNum;
			}
			else
			{
				(*this)(i, j) = 0.0;
			}
		}
	}
}

template <class T>
void  matrix2D<T>::makePositiveDefinite()
{
	this->fillRandomlyLowerTriangular();
	matrix2D<T> transposed = transpose(*this);
	(*this) = multiply((*this), transposed);
}

template <class T>
void matrix2D<T>::resize(const unsigned int newRowCount, const unsigned int newColumnCount)
{
	rowCount = newRowCount;
	columnCount = newColumnCount;
	this->data.resize(rowCount * columnCount);
}

template <class T>
std::ostream& operator<< (std::ostream& stream, const matrix2D<T>& matrix)
{
	matrix.print(stream);
	return stream;
}



template <class T>
matrix2D<T> transpose(matrix2D<T>& matrix)
{
	int n = matrix.columnCount;
	matrix2D<T> matrixT(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			matrixT(i, j) = matrix(j, i);
		}
	}
	return matrixT;
}

template <class T>
bool isValidToMultiply(matrix2D<T> & m1, matrix2D<T> & m2)
{
	bool isValid = true;
	if (m1.columnCount != m2.rowCount)
	{
		std::cerr << "Error in isValidToMultiply(). Num of columns in m1 (" << m1.columnCount << ") is not equal to numbers of rows in m2 ("
			<< m2.rowCount << ")" << std::endl;
		isValid = false;
	}
	return isValid;
}

template <class T>
bool isValidToAdd(matrix2D<T> & m1, matrix2D<T> & m2)
{
	bool isValid = true;
	if (m1.rowCount != m2.rowCount || m1.columnCount != m2.columnCount)
	{
		std::cerr << "Error in isValidToAdd(). Dimensions must match. Got dimensions of : (" << m1.rowCount 
			<< "x" << m1.columnCount << ") and (" << m2.rowCount << "x" << m2.columnCount << ")" << std::endl;
		isValid = false;
	}
	return isValid;
}

template <class T>
matrix2D<T> multiply(matrix2D<T> & m1, matrix2D<T> & m2)
{
	int nRows = m1.rowCount;
	int nCols = m2.columnCount;
	int nColsM1 = m1.columnCount;
	matrix2D<T> result(nRows, nCols);
	if (!isValidToMultiply(m1, m2))
	{
		return result;
	}
	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			// Dot product of ith row from m1 and jth column from m2
			for (int k = 0; k < nColsM1; k++)
			{
				result(i, j) += m1(i, k) * m2(k, j);
				// std::cout << "i=" << i << " j=" << j << " k=" << k << std::endl;
			}
		}
	}
	return result;
}

template <class T>
matrix2D<T> add(matrix2D<T> & m1, matrix2D<T> & m2)
{
	int nRows = m1.rowCount;
	int nCols = m1.columnCount;
	matrix2D<T> result(nRows, nCols);
	if (!isValidToAdd(m1, m2))
	{
		return result;
	}
	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			result(i, j) = m1(i, j) + m2(i, j);
		}
	}
	return result;
}


