#pragma once
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <cmath>
#include <map>

template <class T> class SparseMatrix;

template <class T>
std::ostream& operator<<(std::ostream& stream, const SparseMatrix<T>& matrix);

template <class T>
SparseMatrix<T> transpose(SparseMatrix<T>& matrix);

template <class T>
bool isValidToMultiply(SparseMatrix<T> & m1, SparseMatrix<T> & m2);

template <class T>
bool isValidToAdd(SparseMatrix<T> & m1, SparseMatrix<T> & m2);

//template <class T>
//SparseMatrix<T> multiply(SparseMatrix<T> & m1, SparseMatrix<T> & m2);

//template <class T>
//SparseMatrix<T> add(SparseMatrix<T> & m1, SparseMatrix<T> & m2);

template <class T>
class SparseMatrix
{
public:
	typedef std::multimap<int, std::pair<int, T>> MatrixMap;
	typedef std::pair<int, std::pair<int, T>> MatrixEntry;
	typedef typename MatrixMap::iterator MatrixIter;
	typedef typename MatrixMap::const_iterator MatrixConstIter;
	typedef std::pair<typename MatrixIter, typename MatrixIter> Range;
	SparseMatrix();
	SparseMatrix(int rowCount, int colCount, int nZeroCount=0);
	SparseMatrix(const std::string & mtxFile);
	SparseMatrix(std::vector<T> & v);
	//SparseMatrix(int rowCount, int colCount, T fill);
	//void print() const;
	// int getColumnCount(return columnCount;)
	void print(std::ostream & os) const;
	void insertValue(int row, int col, T val);
	void addToValue(int row, int col, T val);
	// T &operator()(int x, int y);
	// T* getData(int rowIdx);
	// std::vector<T> getVector();
	int getRowCount();
	int getColumnCount();
	int getNumElements();
	int getNumNonZeroElements();
	bool getVal(int row, int col, T & val);
	T & operator()(int row, int col);
	//SparseMatrix<T> createSubMatrix(int startRow, int endRow, int startCol, int endCol);
	//SparseMatrix<T> & getSubMatrix(int startRow, int endRow, int startCol, int endCol);
	//void copyValues(SparseMatrix<T> & from, int startRow, int endRow);
	//void copyValues(SparseMatrix<T> & from, int startRow, int endRow, int startCol, int endCol);
	void fillRandomly(const int min = 0, const int max = 5);
	//void fillRandomlyLowerTriangular(const int min = 0, const int max = 5);
	//void makePositiveDefinite();
	std::vector<T> toVector();
	void fromVector(std::vector<T> & v);
	std::vector<std::vector<T>> getUniqueRows();
	SparseMatrix<T> multiply(SparseMatrix<T> & m1);
	SparseMatrix<T> add(SparseMatrix<T> & m1);
	friend std::ostream& operator<< <>(std::ostream& stream, const SparseMatrix<T>& matrix);
	friend SparseMatrix<T> transpose <>(SparseMatrix<T>& matrix);
	friend bool isValidToMultiply <>(SparseMatrix<T> & m1, SparseMatrix<T> & m2);
	friend bool isValidToAdd <>(SparseMatrix<T> & m1, SparseMatrix<T> & m2);
	typename MatrixMap data;
	//friend SparseMatrix<T> multiply <>(SparseMatrix<T> & m1, SparseMatrix<T> & m2);
	//friend SparseMatrix<T> add <>(SparseMatrix<T> & m1, SparseMatrix<T> & m2);
private:
	//std::vector<T> data;
	int columnCount;
	int rowCount;
	int nonZeroCount;
	//void resize(const unsigned int newRowCount, const unsigned int newColumnCount);
};
template <class T>
SparseMatrix<T>::SparseMatrix() : rowCount(0), columnCount(0), nonZeroCount(0), data()
{
	
}

template <class T>
SparseMatrix<T>::SparseMatrix(int x, int y, int nZeroCount) : rowCount(x), columnCount(y), nonZeroCount(nZeroCount), data()
{

}

template <class T>
SparseMatrix<T>::SparseMatrix(const std::string & mtxFile)
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
			iss >> this->rowCount;
			iss >> this->columnCount;
			//iss >> this->nonZeroCount;
			//this->resize(rowCount, columnCount);
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
			this->insertValue(row, col, val);
		}
	}
	/*
	if (file.bad())
	{
		std::cout << "IO error" << std::endl;// IO error
	}
	else if (!file.eof())
	{
		std::cout << "format error" << std::endl;// format error (not possible with getline but possible with operator>>)
	}
	else 
	{
		std::cout << "format or EOF" << std::endl;
		// format error (not possible with getline but possible with operator>>)
		// or end of file (can't make the difference)
	}
	*/
}
template <class T>
SparseMatrix<T>::SparseMatrix(std::vector<T> & v) : rowCount(0), columnCount(0), nonZeroCount(0), data()
{
	if (v.size() == 0 || v.size() % 3 != 0)
	{
		std::cout << "Warning...created SparseMatrix from an empty or invalid vector" << std::endl;
		return;
	}
	rowCount = (int)v[0];
	columnCount = (int)v[1];
	nonZeroCount = (int)v[2];
	for (int i = 3; i < v.size(); i+=3)
	{
		data.insert(MatrixEntry((int)v[i], std::make_pair((int)v[i + 1], (T)v[i + 2])));
	}
}
/*
template <class T>
SparseMatrix<T>::SparseMatrix(int x, int y, T fill) : data(x*y, fill), rowCount(x), columnCount(y)
{

}
*/
/*
template <class T>
void SparseMatrix<T>::print() const
{
	for (int i = 0; i < rowCount; i++)
	{
		for (int j = 0; j < columnCount; j++)
		{
			//std::cout << data[i * columnCount + j] << ", ";
		}
		std::cout << std::endl;
	}
}
*/
template <class T>
void SparseMatrix<T>::print(std::ostream & os) const
{
	/*
	for (int i = 0; i < rowCount; i++)
	{
		for (int j = 0; j < columnCount; j++)
		{
			//os << data[i * columnCount + j] << ", ";
		}
		os << std::endl;
	}
	*/
	for (MatrixConstIter it = data.begin(); it != data.end(); it++)
	{
		os << "Row: " << it->first << " Column: " << it->second.first 
			<< " Value: " << it->second.second << std::endl;
	}
}

template <class T>
void SparseMatrix<T>::insertValue(int row, int col, T val)
{
	if (val == 0)
	{
		return;
	}
	T oldVal;
	if (this->getVal(row, col, oldVal))
	{
		//oldVal = val;
	}
	else
	{
		// row can col does not exist yet, add it
		data.insert(MatrixEntry(row, std::make_pair(col, val)));
		nonZeroCount++;
	}	
}

template <class T>
void SparseMatrix<T>::addToValue(int row, int col, T val)
{
	if (val == 0)
	{
		return;
	}
	Range r = data.equal_range(row);
	for (MatrixIter lowIt = r.first; lowIt != r.second; lowIt++)
	{
		if (lowIt->second.first == col)
		{
			lowIt->second.second += val;
			return;
		}
	}
	// row can col does not exist yet, add it
	data.insert(MatrixEntry(row, std::make_pair(col, val)));
	nonZeroCount++;
}

/*
template <class T>
T* SparseMatrix<T>::getData(int rowIdx)
{
	return &data[rowIdx * columnCount];
}
*/
/*
template <class T>
std::vector<T> SparseMatrix<T>::getVector()
{
	return data;
}
*/
template <class T>
int SparseMatrix<T>::getRowCount()
{
	return rowCount;
}

template <class T>
int SparseMatrix<T>::getColumnCount()
{
	return columnCount;
}

template <class T>
int SparseMatrix<T>::getNumElements()
{
	return rowCount * columnCount;
}

template <class T>
int SparseMatrix<T>::getNumNonZeroElements()
{
	return nonZeroCount;
}

template <class T>
bool SparseMatrix<T>::getVal(int row, int col, T & val)
{
	if (val == 0)
	{
		return false;
	}
	Range r = data.equal_range(row);
	for (MatrixIter lowIt = r.first; lowIt != r.second; lowIt++)
	{
		if (lowIt->second.first == col)
		{
			val = lowIt->second.second;
			return true;
		}
	}
	return false;
}

template <class T>
T & SparseMatrix<T>::operator()(int row, int col)
{
	Range r = data.equal_range(row);
	for (MatrixIter lowIt = r.first; lowIt != r.second; lowIt++)
	{
		if (lowIt->second.first == col)
		{
			return lowIt->second.second;
		}
	}
	return data.begin()->second.second;
}
/*
template <class T>
SparseMatrix<T> SparseMatrix<T>::createSubMatrix(int startRow, int endRow, int startCol, int endCol)
{
	if (startRow < 0 || startRow > rowCount - 1 || startCol < 0 || startCol > columnCount - 1
		|| endRow < 0 || endRow > rowCount - 1 || endCol < 0 || endCol > columnCount - 1
		|| startRow > endRow || startCol > endCol)
	{
		std::cerr << "createSubMatrix(). Invalid input. startRow=" << startRow << " endRow=" << endRow
			<< " startCol=" << startCol << " endCol=" << endCol << " this->rowCount=" << rowCount << " this->columnCount=" << columnCount << std::endl;
		return SparseMatrix<T>();
	}
	SparseMatrix<T> sub((endRow - startRow) + 1, (endCol - startCol) + 1);
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
*/
/*
template <class T>
SparseMatrix<T> & SparseMatrix<T>::getSubMatrix(int startRow, int endRow, int startCol, int endCol)
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
*/
/*
template <class T>
void  SparseMatrix<T>::copyValues(SparseMatrix<T> & from, int startRow, int endRow)
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
*/
/*
template <class T>
void SparseMatrix<T>::copyValues(SparseMatrix<T> & from, int startRow, int endRow, int startCol, int endCol)
{
	// int numRowsToFillIn = (endRow - startRow) + 1;
	// int numColsToFillIn = (endCol - startCol) + 1;
	int fromRowCount = from.getRowCount();
	int fromColCount = from.getColumnCount();
	// if from is larger than "this" matrix then will fill in as many values as possible from top right
	if (startRow < 0 || startRow > rowCount - 1 || startCol < 0 || startCol > columnCount - 1
		|| endRow < 0 || endRow > rowCount - 1 || endCol < 0 || endCol > columnCount - 1
		|| startRow > endRow || startCol > endCol || from.getColumnCount() > columnCount
												  || from.getRowCount() > rowCount )
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
*/

template <class T>
void SparseMatrix<T>::fillRandomly(const int min, const int max)
{
	/*
	for (int i = 0; i < rowCount; i++)
	{
		for (int j = 0; j < columnCount; j++)
		{
			T randNum = (T)(min + (rand() % static_cast<int>(max - min + 1)));
			randNum = (randNum == 0.0) ? 1.0 : randNum;
			(*this)(i, j) = randNum;
		}
	}
	*/
}
template <class T>
std::vector<T> SparseMatrix<T>::toVector()
{
	std::vector<T> ret;
	ret.push_back(rowCount);
	ret.push_back(columnCount);
	ret.push_back(nonZeroCount);
	for (MatrixConstIter it = data.begin(); it != data.end(); it++)
	{
		ret.push_back(it->first);
		ret.push_back(it->second.first);
		ret.push_back(it->second.second);
	}
	return ret;
}

template <class T>
void SparseMatrix<T>::fromVector(std::vector<T> & v)
{
	if (v.size() == 0 || v.size() % 3 != 0)
	{
		std::cout << "Warning...fromVector() with an empty or invalid vector" << std::endl;
		return;
	}
	data.clear();
	rowCount = (int)v[0];
	columnCount = (int)v[1];
	nonZeroCount = (int)v[2];
	for (int i = 3; i < v.size(); i += 3)
	{
		data.insert(MatrixEntry((int)v[i], std::make_pair((int)v[i + 1], (T)v[i + 2])));
	}
}

template <class T>
std::vector<std::vector<T>> SparseMatrix<T>::getUniqueRows()
{
	std::vector<std::vector<T>> ret;
	for (MatrixConstIter itRows = data.begin(); itRows != data.end(); itRows = data.upper_bound(itRows->first))
	{
		Range r = data.equal_range(itRows->first);
		std::vector<T> currentRow;
		for (MatrixIter it = r.first; it != r.second; it++)
		{
			currentRow.push_back(it->first);
			currentRow.push_back(it->second.first);
			currentRow.push_back(it->second.second);
		}
		ret.push_back(currentRow);
	}
	return ret;
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::multiply(SparseMatrix<T> & m)
{
	int nCols = m.getColumnCount();
	SparseMatrix<T> result(this->rowCount, nCols);
	if (!isValidToMultiply(*this, m))
	{
		return result;
	}
	for (MatrixConstIter it = data.begin(); it != data.end(); it++)
	{
		T val;
		for (int i = 0; i < nCols; i++)
		{
			if (m.getVal(it->second.first, i, val))
			{
				result.addToValue(it->first, i, val * it->second.second);
			}
		}
	}
	return result;
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::add(SparseMatrix<T> & m)
{
	SparseMatrix<T> result(rowCount, columnCount);
	if (!isValidToAdd(*this, m))
	{
		return result;
	}

	for (MatrixConstIter it = data.begin(); it != data.end(); it++)
	{
		T val;
		if (m.getVal(it->first, it->second.first, val))
		{
			result.addToValue(it->first, it->second.first, val + it->second.second);
		}
		else
		{
			result.addToValue(it->first, it->second.first, it->second.second);
		}
	}

	for (MatrixConstIter it = m.data.begin(); it != m.data.end(); it++)
	{
		T val;
		if (!getVal(it->first, it->second.first, val))
		{
			result.addToValue(it->first, it->second.first, it->second.second);
		}
	}
	return result;
}
/*
template <class T>
void SparseMatrix<T>::fillRandomlyLowerTriangular(const int min, const int max)
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
*/
/*
template <class T>
void  SparseMatrix<T>::makePositiveDefinite()
{
	this->fillRandomlyLowerTriangular();
	SparseMatrix<T> transposed = transpose(*this);
	(*this) = multiply((*this), transposed);
}
*/
/*
template <class T>
void SparseMatrix<T>::resize(const unsigned int newRowCount, const unsigned int newColumnCount)
{
	rowCount = newRowCount;
	columnCount = newColumnCount;
	this->data.resize(rowCount * columnCount);
}
*/

template <class T>
std::ostream& operator<< (std::ostream& stream, const SparseMatrix<T>& matrix)
{
	matrix.print(stream);
	return stream;
}

template <class T>
SparseMatrix<T> transpose(SparseMatrix<T>& matrix)
{
	int n = matrix.columnCount;
	SparseMatrix<T> matrixT(n, n);
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
bool isValidToMultiply(SparseMatrix<T> & m1, SparseMatrix<T> & m2)
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
bool isValidToAdd(SparseMatrix<T> & m1, SparseMatrix<T> & m2)
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
/*
template <class T>
SparseMatrix<T> multiply(SparseMatrix<T> & m1, SparseMatrix<T> & m2)
{
	int nRows = m1.rowCount;
	int nCols = m2.columnCount;
	int nColsM1 = m1.columnCount;
	SparseMatrix<T> result(nRows, nCols);
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
	
	//for (MatrixConstIter it = data.begin(); it != data.end(); it++)
	//{

	//}
	return result;
}
*/
/*
template <class T>
SparseMatrix<T> add(SparseMatrix<T> & m1, SparseMatrix<T> & m2)
{
	int nRows = m1.rowCount;
	int nCols = m1.columnCount;
	SparseMatrix<T> result(nRows, nCols);
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
*/
