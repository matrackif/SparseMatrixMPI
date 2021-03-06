﻿#pragma once
#include "matrix.hpp"
#include "mpi.h"

template <typename T>
void concat(std::vector<T>& a, const std::vector<T>& b)
{
	a.reserve(a.size() + b.size());
	a.insert(a.end(), b.begin(), b.end());
}
template <typename T>
int getNumElementsInRange(const std::vector<T> & v, const int & startRow, const int & endRow, int & startIdx)
{
	int ret = 0;
	int vecSize = v.size();
	bool assignedStartIdx = false;
	for (int i = 3; i < vecSize; i += 3)
	{
		if (v[i] >= startRow && v[i] < endRow)
		{
			if (!assignedStartIdx)
			{
				startIdx = i;
				assignedStartIdx = true;
			}
			ret += 3;
		}
		else if (v[i] >= endRow)
		{
			break;
		}

	}
	return ret;
}
// TODO use block striped decomposition algo http://www.hpcc.unn.ru/mskurs/ENG/PPT/pp08.pdf
template <class T>
void parallelMult(int rank, int size, matrix2D<T> & m1, matrix2D<T> & m2)
{
	const int tag = 0;
	MPI_Status status;
	// TODO have rank 0 bcast matrix dims to all other ranks
	int m1RowCount = m1.getRowCount();
	int m1ColCount = m1.getColumnCount();
	int m2RowCount = m2.getRowCount();
	int m2ColCount = m2.getColumnCount();
	MPI_Bcast(&m1RowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m1ColCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m2RowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m2ColCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// std::cout << "Rank=" << rank << " m1RowCount: " << m1RowCount << " m1ColCount: " << m1ColCount << " m2RowCount: " << m2RowCount << " m2ColCount: " << m2ColCount << std::endl;
	// std::cout << "Rank=" << rank << " matrix m1: \n" << m1 << std::endl;
	// std::cout << "Rank=" << rank << " matrix m2: \n" << m2 << std::endl;
	if (m1RowCount <= 0 || m1ColCount <= 0 || m2RowCount <= 0 || m2ColCount <= 0)
	{
		if (rank == 0)
		{
			std::cerr << "parallelMult() Error: One of the matrix dimensions is <= 0" << std::endl;
		}
		return;
	}


	// matrix2D<T> result, m1SubMatrix, m2SubMatrix;
	matrix2D<T> result, m1SubMatrix, m2SubMatrix;
	int m1MinRowsPerProcess = m1RowCount / size;
	int m1RowRemainder = m1RowCount % size;
	int numRowsPerProcess = rank < (m1RowRemainder) ? m1MinRowsPerProcess + 1 : m1MinRowsPerProcess;
	if (rank == 0)
	{
		// int offset = (m1RowCount % size) + (m1RowCount / size);
		//std::cout << "m1RowCount: " << m1RowCount << " m1ColCount: " << m1ColCount << " m2RowCount: " << m2RowCount << " m2ColCount: " << m2ColCount << std::endl;
		std::cout << "parallelMult() Rank=" << rank << " numRowsPerProcess: " << numRowsPerProcess << std::endl;
		//std::cout << "matrix m1: \n" << m1 << std::endl;
		//std::cout << "matrix m2: \n" << m2 << std::endl;
		result = matrix2D<T>(m1RowCount, m2ColCount);
		// Send rest of matrix, then calculate
		int numRowsToSend;
		for (int d = 0; d < size - 1; d++)
		{
			numRowsToSend = (d + 1) < m1RowRemainder ? m1MinRowsPerProcess + 1 : m1MinRowsPerProcess; // we have to recalculate amount of rows a specific rank has?
			MPI_Send(m1.getData(numRowsPerProcess + (d * numRowsToSend)), m1ColCount * numRowsToSend, MPI_DOUBLE, d + 1, tag, MPI_COMM_WORLD);
			MPI_Send(m2.getData(0), m2.getNumElements(), MPI_DOUBLE, d + 1, tag, MPI_COMM_WORLD); // TODO partition it differently
		}
		// matrix2D<T> & m1SubMatrix = m1.getSubMatrix(0, numRowsPerProcess - 1, 0, m1ColCount - 1) ALTERNATE METHOD THAT MODIFIES m1 but is faster!!
		m1SubMatrix = m1.createSubMatrix(0, numRowsPerProcess - 1, 0, m1ColCount - 1); // NOTE: using getSubMatrix()
		m2SubMatrix = m2.createSubMatrix(0, m2RowCount - 1, 0, m2ColCount - 1); // TODO partition it differently
		// std::cout << "Rank="<< rank << " m1 submatrix: \n" << m1SubMatrix << std::endl;
		// std::cout << "Rank=" << rank << " m2 submatrix: \n" << m2SubMatrix << std::endl;
		matrix2D<T> partialResult = multiply(m1SubMatrix, m2SubMatrix);
		// std::cout << "Rank=" << rank << " multiply() m1SubMatrix* m2SubMatrix \n" << partialResult << std::endl;
		result.copyValues(partialResult, 0, numRowsPerProcess - 1);
		// std::cout << "Rank=" << rank << " copyValues() done" << std::endl;
		int numRowsToReceive;
		for (int i = 0; i < size - 1; i++)
		{
			numRowsToReceive = (i + 1) < m1RowRemainder ? m1MinRowsPerProcess + 1 : m1MinRowsPerProcess; // we have to recalculate amount of rows a specific rank has?
			MPI_Recv(result.getData(numRowsPerProcess + (i * numRowsToReceive)), numRowsToReceive * m2ColCount, MPI_DOUBLE, i + 1, tag, MPI_COMM_WORLD, &status);
		}
		std::cout << "parallelMult() Finished! m1: \n" << m1 << "m2: \n" << m2 << " m1 * m2: \n" << result << std::endl;
	}
	else
	{
		std::cout << "parallelMult() Rank=" << rank << " numRowsPerProcess: " << numRowsPerProcess << std::endl;
		result = matrix2D<T>(numRowsPerProcess, m2ColCount);
		m1SubMatrix = matrix2D<T>(numRowsPerProcess, m1ColCount);
		m2SubMatrix = matrix2D<T>(m2RowCount, m2ColCount); // TODO partition it differently
		MPI_Recv(m1SubMatrix.getData(0), m1SubMatrix.getNumElements(), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(m2SubMatrix.getData(0), m2SubMatrix.getNumElements(), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		// std::cout << "Rank=" << rank << " m1 submatrix: \n" << m1SubMatrix << std::endl;
		// std::cout << "Rank=" << rank << " m2 submatrix: \n" << m2SubMatrix << std::endl;
		matrix2D<T> partialResult = multiply(m1SubMatrix, m2SubMatrix);
		MPI_Send(partialResult.getData(0), partialResult.getNumElements(), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}
}

template <class T>
void parallelAdd(int rank, int size, matrix2D<T> & m1, matrix2D<T> & m2)
{
	const int tag = 0;
	MPI_Status status;
	// TODO have rank 0 bcast matrix dims to all other ranks
	int m1RowCount = m1.getRowCount();
	int m1ColCount = m1.getColumnCount();
	int m2RowCount = m2.getRowCount();
	int m2ColCount = m2.getColumnCount();
	MPI_Bcast(&m1RowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m1ColCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m2RowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m2ColCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// std::cout << "Rank=" << rank << " m1RowCount: " << m1RowCount << " m1ColCount: " << m1ColCount << " m2RowCount: " << m2RowCount << " m2ColCount: " << m2ColCount << std::endl;
	// std::cout << "Rank=" << rank << " matrix m1: \n" << m1 << std::endl;
	// std::cout << "Rank=" << rank << " matrix m2: \n" << m2 << std::endl;
	if (m1RowCount != m2RowCount || m1ColCount != m2ColCount)
	{
		if (rank == 0)
		{
			std::cerr << "parallelAdd() Error: Dimensions must be equal" << std::endl;
		}
		return;
	}

	matrix2D<T> result, m1SubMatrix, m2SubMatrix;
	int m1MinRowsPerProcess = m1RowCount / size;
	int m1RowRemainder = m1RowCount % size;
	int numRowsPerProcess = rank < (m1RowRemainder) ? m1MinRowsPerProcess + 1 : m1MinRowsPerProcess;
	if (rank == 0)
	{
		// int offset = (m1RowCount % size) + (m1RowCount / size);
		//std::cout << "m1RowCount: " << m1RowCount << " m1ColCount: " << m1ColCount << " m2RowCount: " << m2RowCount << " m2ColCount: " << m2ColCount << std::endl;
		std::cout << "parallelAdd() Rank=" << rank << " numRowsPerProcess: " << numRowsPerProcess << std::endl;
		//std::cout << "matrix m1: \n" << m1 << std::endl;
		//std::cout << "matrix m2: \n" << m2 << std::endl;
		result = matrix2D<T>(m1RowCount, m1ColCount);
		// Send rest of matrix, then calculate
		int numRowsToSend;
		for (int d = 0; d < size - 1; d++)
		{
			numRowsToSend = (d + 1) < m1RowRemainder ? m1MinRowsPerProcess + 1 : m1MinRowsPerProcess; 
			MPI_Send(m1.getData(numRowsPerProcess + (d * numRowsToSend)), m1ColCount * numRowsToSend, MPI_DOUBLE, d + 1, tag, MPI_COMM_WORLD);
			MPI_Send(m2.getData(numRowsPerProcess + (d * numRowsToSend)), m2ColCount * numRowsToSend, MPI_DOUBLE, d + 1, tag, MPI_COMM_WORLD); // TODO partition it differently
		}
		m1SubMatrix = m1.createSubMatrix(0, numRowsPerProcess - 1, 0, m1ColCount - 1);
		m2SubMatrix = m2.createSubMatrix(0, numRowsPerProcess - 1, 0, m2ColCount - 1);
		// std::cout << "Rank=" << rank << " m1 submatrix: \n" << m1SubMatrix << std::endl;
		// std::cout << "Rank=" << rank << " m2 submatrix: \n" << m2SubMatrix << std::endl;
		matrix2D<T> partialResult = add(m1SubMatrix, m2SubMatrix);
		// std::cout << "Rank=" << rank << " add() m1SubMatrix + m2SubMatrix \n" << partialResult << std::endl;
		result.copyValues(partialResult, 0, numRowsPerProcess - 1);
		// std::cout << "Rank=" << rank << " copyValues() done" << std::endl;
		int numRowsToReceive;
		for (int i = 0; i < size - 1; i++)
		{
			numRowsToReceive = (i + 1) < m1RowRemainder ? m1MinRowsPerProcess + 1 : m1MinRowsPerProcess;
			MPI_Recv(result.getData(numRowsPerProcess + (i * numRowsToReceive)), numRowsToReceive * m2ColCount, MPI_DOUBLE, i + 1, tag, MPI_COMM_WORLD, &status);
		}
		std::cout << "parallelAdd() Finished! m1: \n" << m1 << "m2: \n" << m2 << " m1 + m2: \n" << result << std::endl;
	}
	else 
	{
		std::cout << "parallelAdd() Rank=" << rank << " numRowsPerProcess: " << numRowsPerProcess << std::endl;
		result = matrix2D<T>(numRowsPerProcess, m1ColCount);
		m1SubMatrix = matrix2D<T>(numRowsPerProcess, m1ColCount);
		m2SubMatrix = matrix2D<T>(numRowsPerProcess, m2ColCount); // TODO partition it differently
		MPI_Recv(m1SubMatrix.getData(0), m1SubMatrix.getNumElements(), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(m2SubMatrix.getData(0), m2SubMatrix.getNumElements(), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		// std::cout << "Rank=" << rank << " m1 submatrix: \n" << m1SubMatrix << std::endl;
		// std::cout << "Rank=" << rank << " m2 submatrix: \n" << m2SubMatrix << std::endl;
		matrix2D<T> partialResult = add(m1SubMatrix, m2SubMatrix);
		MPI_Send(partialResult.getData(0), partialResult.getNumElements(), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}
}


template <class T>
SparseMatrix<T> parallelAdd(int rank, int size, SparseMatrix<T> & m1, SparseMatrix<T> & m2)
{
	const int tag = 0;
	MPI_Status status;
	// TODO have rank 0 bcast matrix dims to all other ranks
	std::vector<std::vector<double> > uniqueRows;
	std::vector<std::vector<double> > uniqueRows2;
	if (rank == 0)
	{
		uniqueRows = m1.getUniqueRows();
		uniqueRows2 = m2.getUniqueRows();
	}
	int m1UniqueRowCount = uniqueRows.size();
	int m2UniqueRowCount = uniqueRows2.size();
	int m1RowCount = m1.getRowCount();
	int m1ColCount = m1.getColumnCount();
	int m2RowCount = m2.getRowCount();
	int m2ColCount = m2.getColumnCount();
	MPI_Bcast(&m1UniqueRowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m2UniqueRowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m1RowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m1ColCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m2RowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m2ColCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// std::cout << "Rank=" << rank << " m1RowCount: " << m1RowCount << " m1ColCount: " << m1ColCount << " m2RowCount: " << m2RowCount << " m2ColCount: " << m2ColCount << std::endl;
	// std::cout << "Rank=" << rank << " matrix m1: \n" << m1 << std::endl;
	// std::cout << "Rank=" << rank << " matrix m2: \n" << m2 << std::endl;
	// matrix2D<T> result, m1SubMatrix, m2SubMatrix;
	if (m1RowCount != m2RowCount || m1ColCount != m2ColCount)
	{
		if (rank == 0)
		{
			std::cerr << "Sparse Matrix parallelAdd(). Error, dimensions don't match" << std::endl;
		}
		return SparseMatrix<T>();
	}
	SparseMatrix<T> result, m1SubMatrix, m2SubMatrix;
	int m1MinRowsPerProcess = m1UniqueRowCount / size;
	int m1RowRemainder = m1UniqueRowCount % size;
	int numRowsPerProcess = rank < (m1RowRemainder) ? m1MinRowsPerProcess + 1 : m1MinRowsPerProcess;

	int m2MinRowsPerProcess = m2UniqueRowCount / size;
	int m2RowRemainder = m2UniqueRowCount % size;
	int numRowsPerProcessM2 = rank < (m2RowRemainder) ? m2MinRowsPerProcess + 1 : m2MinRowsPerProcess;
	if (rank == 0)
	{
		// int offset = (m1RowCount % size) + (m1RowCount / size);
		//std::cout << "m1RowCount: " << m1RowCount << " m1ColCount: " << m1ColCount << " m2RowCount: " << m2RowCount << " m2ColCount: " << m2ColCount << std::endl;
		//std::cout << "parallelAdd() Rank=" << rank << " numRowsPerProcess: " << numRowsPerProcess << std::endl;
		//std::cout << "parallelAdd() Rank=" << rank << " numRowsPerProcessM2: " << numRowsPerProcessM2 << std::endl;
		//std::cout << "matrix m1: \n" << m1 << std::endl;
		//std::cout << "matrix m2: \n" << m2 << std::endl;
		result = SparseMatrix<T>(m1.getRowCount(), m2.getColumnCount());
		std::vector<T> resultAsVec = result.toVector();
		//std::cout << "resultAsVec: \n";
		/*
		for (int k = 0; k < resultAsVec.size(); k++)
		{
			std::cout << resultAsVec[k] << ",";
		}
		*/
		std::vector<T> m2AsVec = m2.toVector();
		// Send rest of matrix, then calculate
		int numRowsToSend, numRowsToSendM2, uniqueRowIdx = numRowsPerProcess, uniqueRowIdx2 = numRowsPerProcessM2;
		for (int d = 0; d < size - 1; d++)
		{
			numRowsToSend = (d + 1) < m1RowRemainder ? m1MinRowsPerProcess + 1 : m1MinRowsPerProcess; // we have to recalculate amount of rows a specific rank has?
			numRowsToSendM2 = (d + 1) < m2RowRemainder ? m2MinRowsPerProcess + 1 : m2MinRowsPerProcess;
			std::vector<T> currentMatrixVec;
			currentMatrixVec.push_back(m1.getRowCount()); 
			currentMatrixVec.push_back(m1ColCount);
			currentMatrixVec.push_back(0); 

			std::vector<T> currentMatrixVec2;
			currentMatrixVec2.push_back(m1.getRowCount());
			currentMatrixVec2.push_back(m1ColCount);
			currentMatrixVec2.push_back(0);
			for (int i = 0; i < numRowsToSend; i++)
			{
				concat(currentMatrixVec, uniqueRows[uniqueRowIdx]);
				++uniqueRowIdx;
			}
			//std::cout << "numRowsToSendM2: " << numRowsToSendM2 <<std::endl;
			for (int i = 0; i < numRowsToSendM2; i++)
			{
				concat(currentMatrixVec2, uniqueRows2[uniqueRowIdx2]);
				++uniqueRowIdx2;
			}
			int m1Size = currentMatrixVec.size();
			int m2Size = currentMatrixVec2.size();
			MPI_Send(&m1Size, 1, MPI_INT, d + 1, tag, MPI_COMM_WORLD);
			MPI_Send(&m2Size, 1, MPI_INT, d + 1, tag, MPI_COMM_WORLD);
			MPI_Send(&currentMatrixVec[0], m1Size, MPI_DOUBLE, d + 1, tag, MPI_COMM_WORLD);
			MPI_Send(&currentMatrixVec2[0], m2Size, MPI_DOUBLE, d + 1, tag, MPI_COMM_WORLD); // TODO partition it differently
		}
		std::vector<T> m1SubMatrixVec;
		m1SubMatrixVec.push_back(m1RowCount); // dummy value
		m1SubMatrixVec.push_back(m1ColCount);
		m1SubMatrixVec.push_back(0); // dummy value

		std::vector<T> m2SubMatrixVec;
		m2SubMatrixVec.push_back(m1.getRowCount()); // dummy value
		m2SubMatrixVec.push_back(m1ColCount);
		m2SubMatrixVec.push_back(0); // dummy value
		for (int i = 0; i < numRowsPerProcess; i++)
		{
			concat(m1SubMatrixVec, uniqueRows[i]);
		}
		for (int i = 0; i < numRowsPerProcessM2; i++)
		{
			concat(m2SubMatrixVec, uniqueRows2[i]);
		}
		/*
		std::cout << "m1SubMatrixVec: \n";
		for (int k = 0; k < m1SubMatrixVec.size(); k++)
		{
			std::cout << m1SubMatrixVec[k] << ",";
		}
		*/
		m1SubMatrix.fromVector(m1SubMatrixVec);
		m2SubMatrix.fromVector(m2SubMatrixVec);
		//std::cout << "m1SubMatrix: \n" << m1SubMatrix << std::endl;
		SparseMatrix<T> partialResult = m1SubMatrix.add(m2SubMatrix);
		//std::cout << "m2: \n" << m2 << std::endl;
		//std::cout << "partialResult: \n" << partialResult << std::endl;
		std::vector<T> partialResultVec = partialResult.toVector();
		for (int i = 3; i < partialResultVec.size(); i++)
		{
			resultAsVec.push_back(partialResultVec[i]);
		}
		//int numRowsToReceive;
		for (int i = 0; i < size - 1; i++)
		{
			//numRowsToReceive = (i + 1) < m1RowRemainder ? m1MinRowsPerProcess + 1 : m1MinRowsPerProcess; // we have to recalculate amount of rows a specific rank has?
			int size;
			MPI_Recv(&size, 1, MPI_INT, i + 1, tag, MPI_COMM_WORLD, &status);
			//std::cout << "Rank " << i + 1 << " has sent size: " << size << std::endl;
			std::vector<T> tmp(size);
			MPI_Recv(&tmp[0], size, MPI_DOUBLE, i + 1, tag, MPI_COMM_WORLD, &status);
			for (int j = 3; j < size; j++)
			{
				//std::cout << "Pushing" << tmp[j] << std::endl;
				resultAsVec.push_back(tmp[j]);
			}
		}
		/*
		std::cout << "resultAsVec: \n";
		for (int k = 0; k < resultAsVec.size(); k++)
		{
			std::cout << resultAsVec[k] << ",";
		}
		*/
		result.fromVector(resultAsVec);
		//std::cout << "SparseMatrix parallelAdd() Finished! m1: \n" << m1 << "m2: \n" << m2 << " m1 + m2: \n" << result << std::endl;
		return result;
	}
	else
	{
		//std::cout << "parallelAdd() Rank=" << rank << " numRowsPerProcess: " << numRowsPerProcess << std::endl;
		//std::cout << "parallelAdd() Rank=" << rank << " numRowsPerProcessM2: " << numRowsPerProcessM2 << std::endl;
		int m1Size, m2Size;
		MPI_Recv(&m1Size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&m2Size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		//std::cout << "parallelAdd() Rank=" << rank << " m1Size: " << m1Size << " m2Size: " << m2Size << std::endl;
		std::vector<T> m1SubMatrixVec(m1Size), m2SubMatrixVec(m2Size);
		MPI_Recv(&m1SubMatrixVec[0], m1Size, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&m2SubMatrixVec[0], m2Size, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		// std::cout << "Rank=" << rank << " m1 submatrix: \n" << m1SubMatrix << std::endl;
		// std::cout << "Rank=" << rank << " m2 submatrix: \n" << m2SubMatrix << std::endl;
		m1SubMatrix = SparseMatrix<T>(m1SubMatrixVec);
        m2SubMatrix = SparseMatrix<T>(m2SubMatrixVec);
		std::vector<T> partialResultVec = (m1SubMatrix.add(m2SubMatrix)).toVector();
		int size = partialResultVec.size();
		MPI_Send(&size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
		MPI_Send(&partialResultVec[0], size, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}
	return SparseMatrix<T>();
}

template <class T>
SparseMatrix<T> parallelMult(int rank, int size, SparseMatrix<T> & m1, SparseMatrix<T> & m2)
{
	const int tag = 0;
	MPI_Status status;
	// TODO have rank 0 bcast matrix dims to all other ranks
	std::vector<std::vector<double> > uniqueRows;
	if (rank == 0)
	{
		uniqueRows = m1.getUniqueRows();
	}
	int m1UniqueRowCount = uniqueRows.size();
	int m1ColCount = m1.getColumnCount();
	int m2RowCount = m2.getRowCount();
	int m2ColCount = m2.getColumnCount();
	MPI_Bcast(&m1UniqueRowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m1ColCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m2RowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m2ColCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// std::cout << "Rank=" << rank << " m1RowCount: " << m1RowCount << " m1ColCount: " << m1ColCount << " m2RowCount: " << m2RowCount << " m2ColCount: " << m2ColCount << std::endl;
	// std::cout << "Rank=" << rank << " matrix m1: \n" << m1 << std::endl;
	// std::cout << "Rank=" << rank << " matrix m2: \n" << m2 << std::endl;
	/*
	if (m1RowCount <= 0 || m1ColCount <= 0 || m2RowCount <= 0 || m2ColCount <= 0)
	{
	if (rank == 0)
	{
	std::cerr << "parallelMult() Error: One of the matrix dimensions is <= 0" << std::endl;
	}
	return;
	}
	*/

	// matrix2D<T> result, m1SubMatrix, m2SubMatrix;
	SparseMatrix<T> result, m1SubMatrix, m2SubMatrix;
	int m1MinRowsPerProcess = m1UniqueRowCount / size;
	int m1RowRemainder = m1UniqueRowCount % size;
	int numRowsPerProcess = rank < (m1RowRemainder) ? m1MinRowsPerProcess + 1 : m1MinRowsPerProcess;
	if (rank == 0)
	{
		// int offset = (m1RowCount % size) + (m1RowCount / size);
		//std::cout << "m1RowCount: " << m1RowCount << " m1ColCount: " << m1ColCount << " m2RowCount: " << m2RowCount << " m2ColCount: " << m2ColCount << std::endl;
		//std::cout << "parallelMult() Rank=" << rank << " numRowsPerProcess: " << numRowsPerProcess << std::endl;
		//std::cout << "matrix m1: \n" << m1 << std::endl;
		//std::cout << "matrix m2: \n" << m2 << std::endl;
		result = SparseMatrix<T>(m1.getRowCount(), m2.getColumnCount());
		std::vector<T> resultAsVec = result.toVector();
		//std::cout << "resultAsVec: \n";
		//for (int k = 0; k < resultAsVec.size(); k++)
		//{
		//	std::cout << resultAsVec[k] << ",";
		//}
		std::vector<T> m2AsVec = m2.toVector();
		// Send rest of matrix, then calculate
		int numRowsToSend, uniqueRowIdx = numRowsPerProcess;
		for (int d = 0; d < size - 1; d++)
		{
			numRowsToSend = (d + 1) < m1RowRemainder ? m1MinRowsPerProcess + 1 : m1MinRowsPerProcess; // we have to recalculate amount of rows a specific rank has?
			std::vector<T> currentMatrixVec;
			currentMatrixVec.push_back(m1.getRowCount()); // dummy value
			currentMatrixVec.push_back(m1ColCount);
			currentMatrixVec.push_back(0); // dummy value
			for (int i = 0; i < numRowsToSend; i++)
			{
				concat(currentMatrixVec, uniqueRows[uniqueRowIdx]);
				++uniqueRowIdx;
			}
			int m1Size = currentMatrixVec.size();
			int m2Size = m2AsVec.size();
			MPI_Send(&m1Size, 1, MPI_INT, d + 1, tag, MPI_COMM_WORLD);
			MPI_Send(&m2Size, 1, MPI_INT, d + 1, tag, MPI_COMM_WORLD);
			MPI_Send(&currentMatrixVec[0], m1Size, MPI_DOUBLE, d + 1, tag, MPI_COMM_WORLD);
			MPI_Send(&m2AsVec[0], m2Size, MPI_DOUBLE, d + 1, tag, MPI_COMM_WORLD); // TODO partition it differently
		}
		std::vector<T> m1SubMatrixVec;
		m1SubMatrixVec.push_back(m1.getRowCount()); // dummy value
		m1SubMatrixVec.push_back(m1ColCount);
		m1SubMatrixVec.push_back(0); // dummy value
		for (int i = 0; i < numRowsPerProcess; i++)
		{
			concat(m1SubMatrixVec, uniqueRows[i]);
		}
		/*
		std::cout << "m1SubMatrixVec: \n";
		for (int k = 0; k < m1SubMatrixVec.size(); k++)
		{
		std::cout << m1SubMatrixVec[k] << ",";
		}
		*/
		m1SubMatrix.fromVector(m1SubMatrixVec);
		//std::cout << "m1SubMatrix: \n" << m1SubMatrix << std::endl;
		SparseMatrix<T> partialResult = m1SubMatrix.multiply(m2);
		//std::cout << "m2: \n" << m2 << std::endl;
		//std::cout << "partialResult: \n" << partialResult << std::endl;
		std::vector<T> partialResultVec = partialResult.toVector();
		for (int i = 3; i < partialResultVec.size(); i++)
		{
			resultAsVec.push_back(partialResultVec[i]);
		}
		int numRowsToReceive;
		for (int i = 0; i < size - 1; i++)
		{
			numRowsToReceive = (i + 1) < m1RowRemainder ? m1MinRowsPerProcess + 1 : m1MinRowsPerProcess; // we have to recalculate amount of rows a specific rank has?
			int size;
			MPI_Recv(&size, 1, MPI_INT, i + 1, tag, MPI_COMM_WORLD, &status);
			//std::cout << "Rank " << i + 1 << " has sent size: " << size << std::endl;
			std::vector<T> tmp(size);
			MPI_Recv(&tmp[0], size, MPI_DOUBLE, i + 1, tag, MPI_COMM_WORLD, &status);
			for (int j = 3; j < size; j++)
			{
				//std::cout << "Pushing" << tmp[j] << std::endl;
				resultAsVec.push_back(tmp[j]);
			}
		}
		/*
		std::cout << "resultAsVec: \n";
		for (int k = 0; k < resultAsVec.size(); k++)
		{
		std::cout << resultAsVec[k] << ",";
		}
		*/
		result.fromVector(resultAsVec);
		//std::cout << "SparseMatrix parallelMult() Finished! m1: \n" << m1 << "m2: \n" << m2 << " m1 * m2: \n" << result << std::endl;
		return result;
	}
	else
	{
		//std::cout << "parallelMult() Rank=" << rank << " numRowsPerProcess: " << numRowsPerProcess << std::endl;

		int m1Size, m2Size;
		MPI_Recv(&m1Size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&m2Size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		//std::cout << "parallelMult() Rank=" << rank << " m1Size: " << m1Size << " m2Size: " << m2Size << std::endl;
		std::vector<T> m1SubMatrixVec(m1Size), m2SubMatrixVec(m2Size);
		MPI_Recv(&m1SubMatrixVec[0], m1Size, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&m2SubMatrixVec[0], m2Size, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		// std::cout << "Rank=" << rank << " m1 submatrix: \n" << m1SubMatrix << std::endl;
		// std::cout << "Rank=" << rank << " m2 submatrix: \n" << m2SubMatrix << std::endl;
        m1SubMatrix = SparseMatrix<T>(m1SubMatrixVec);
        m2SubMatrix = SparseMatrix<T>(m2SubMatrixVec);
		std::vector<T> partialResultVec = (m1SubMatrix.multiply(m2SubMatrix)).toVector();
		int size = partialResultVec.size();
		MPI_Send(&size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
		MPI_Send(&partialResultVec[0], size, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}
	return SparseMatrix<T>();
}

template <class T>
void conjugateGradient(int rank, int size, SparseMatrix<T> & A, SparseMatrix<T> & b, SparseMatrix<T> & x)
{
	int maxIters = b.getRowCount();
	MPI_Bcast(&maxIters, 1, MPI_INT, 0, MPI_COMM_WORLD);

	SparseMatrix<T> NegativeAx = -1 * parallelMult(rank, size, A, x);
	SparseMatrix<T> r = parallelAdd(rank, size, b, NegativeAx);
	SparseMatrix<T> p = r;
	SparseMatrix<T> rTransposed = transpose(r);
	SparseMatrix<T> rsold = parallelMult(rank, size, rTransposed, r);
	T rsoldT = rsold.getFirstVal();
	/*
	if (rank == 0)
	{
		std::cout << "r and p: \n" << r << std::endl;
		std::cout << "rsold: \n" << rsold << std::endl;
		std::cout << "rsoldT: " << rsoldT << std::endl;
	}
	*/
	for (int i = 0; i < maxIters; i++)
	{
		SparseMatrix<T> Ap = parallelMult(rank, size, A, p);
		SparseMatrix<T> pTransposed = transpose(p);
		T alpha = rsoldT / (parallelMult(rank, size, pTransposed, Ap)).getFirstVal();
		//SparseMatrix<T> alphaSm(1, 1, 0);
		//alphaSm.insertValue(0, 0, rsoldT);
		// TODO maybe just let rank 0 multiply a scalar by a matrix. instead of doing it in parallel?
		//std::cout << "iter: " << i << " Alpha * p: \n" << alpha * p << std::endl;
		SparseMatrix<T> alphaP = alpha * p;
		x = parallelAdd(rank, size, x, alphaP);
		SparseMatrix<T> negativeAlphaAp = -1 * (alpha * Ap);
		r = parallelAdd(rank, size, r, negativeAlphaAp);
		SparseMatrix<T>  newRTransposed = transpose(r);
		SparseMatrix<T> rsnew = parallelMult(rank, size, newRTransposed, r);
		T rsnewT = rsnew.getFirstVal();
		SparseMatrix<T> newP = (rsnewT / rsoldT) * p;
		p = parallelAdd(rank, size, r, newP);
		rsoldT = rsnewT;
		/*
		if (rank == 0)
		{
			std::cout << "iter: " << i << " Ap: \n" << Ap << std::endl;
			std::cout << "iter: " << i << " alpha: \n" << alpha << std::endl;
			std::cout << "iter: " << i << " x: \n" << x << std::endl;
			std::cout << "iter: " << i << " r: \n" << r << std::endl;
			std::cout << "iter: " << i << " rsnew: \n" << rsnew << std::endl;
		}
		*/
	}
	/*
	if (rank == 0)
	{
		std::cout << "conjugateGradient() finished, x is: \n" << x << std::endl;
	}
	*/
}

template <class T>
void incompleteCholeskyDecomp(SparseMatrix<T> & a)
{
	int n = a.getRowCount(); // assume square matrix
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	for (int k = 0; k < n; k++)
	{
		typename std::multimap<int, std::pair<int, T> >::iterator itKK = a.getIter(k, k);
		if (itKK != a.data.end())
		{
			itKK->second.second = sqrt(itKK->second.second);
			for (int i = k + 1; i < n; i++)
			{
				typename std::multimap<int, std::pair<int, T> >::iterator itIK = a.getIter(i, k);
				if (itIK != a.data.end())
				{
					itIK->second.second = itIK->second.second / itKK->second.second;
				}
			}
		}

		for (int j = k + 1; j < n; j++)
		{
			for (int i = j; i < n; i++)
			{
				typename std::multimap<int, std::pair<int, T> >::iterator itIJ = a.getIter(i, j);
				if (itIJ != a.data.end())
				{
					typename std::multimap<int, std::pair<int, T> >::iterator itIK = a.getIter(i, k);
					typename std::multimap<int, std::pair<int, T> >::iterator itJK = a.getIter(j, k);
					if (itIK == a.data.end() || itJK == a.data.end())
					{
						a.data.erase(itIJ);
					}
					else
					{
						itIJ->second.second -= itIK->second.second * itJK->second.second;
					}
				}
			}
		}
		
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			a.deleteVal(i, j);
		}
	}
}
/*
for k = 1, . . . , n do
	for i = k + 1, . . . , n do
		aik = aik /akk
		for j = k + 1, . . . , n do
			aij = aij − aik akj
		end for
	end for
end for
*/

template <class T>
void parallelILU(int rank, int size, SparseMatrix<T> & a)
{
	int n = a.getRowCount(), vecSize = 0; // assume square matrix
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	std::vector<T> asVec;
	SparseMatrix<T> ret;
	if(rank == 0)
	{
		asVec = a.toVector();
		vecSize = asVec.size();
	}
	MPI_Bcast(&vecSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank != 0)
	{
		asVec.resize(vecSize);
	}
	//std::cout << "Rank: " << rank << " sees vecSize as: " << vecSize << "but vec.size() as: " << asVec.size() << std::endl;
	MPI_Bcast(&asVec[0], vecSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	int minRowsPerProcess = n / size;
	int rowRemainder = n % size;
	int numRowsPerProcess = rank < (rowRemainder) ? minRowsPerProcess + 1 : minRowsPerProcess;
	int startRow = (rank >= rowRemainder) ? ((rowRemainder * (minRowsPerProcess + 1)) + ((rank - rowRemainder) *  minRowsPerProcess)) : (numRowsPerProcess * rank);
	int endRow = startRow + numRowsPerProcess;
	//std::cout << "Rank: " << rank << " startRow: " << startRow << " endRow: " << endRow << std::endl;
	/*
	for (auto & lol : asVec)
	{
		std::cout << lol << ", ";
	}
	std::cout << std::endl;
	*/
	
	for (int k = 0; k < n; k++)
	{
		SparseMatrix<T> localCopy(asVec);
		int startIdx = -1;
		int numElements = getNumElementsInRange(asVec, startRow, endRow, startIdx);
		if (startIdx > 0 && numElements > 0)
		{
			startRow = startRow < k + 1 ? k + 1 : startRow;
			for (int i = startRow; i < endRow; i++)
			{
				typename std::multimap<int, std::pair<int, T> >::iterator itKK = localCopy.getIter(k, k);
				//if (itKK != localCopy.data.end())
				//{
				typename std::multimap<int, std::pair<int, T> >::iterator itIK = localCopy.getIter(i, k);
				if (itIK != localCopy.data.end())
				{
					itIK->second.second /= itKK->second.second;
					for (int j = k + 1; j < n; j++)
					{
						typename std::multimap<int, std::pair<int, T> >::iterator itKJ = localCopy.getIter(k, j);
						if (itKJ != localCopy.data.end())
						{
							typename std::multimap<int, std::pair<int, T> >::iterator itIJ = localCopy.getIter(i, j);
							if (itIJ == localCopy.data.end())
							{
								// insertValue checks again
								localCopy.data.insert(std::pair<int, std::pair<int, T> >(i, std::make_pair(j, -1 * (itIK->second.second * itKJ->second.second))));
								numElements += 3;
								std::cout << "adding new pair" << std::endl;
							}
							else 
							{
								itIJ->second.second -= -1 * (itIK->second.second * itKJ->second.second);
								//std::cout << "Rank: " << rank << "modifying ("<< i << "," << j << ")  to " << itIJ->second.second << std::endl;
							}
						}
					}
				}
				//}
			}
			
		}
		asVec = localCopy.toVector();
		if (startIdx < 0)
		{
			//std::cout << "Rank: " << rank << " iter k: " << k << " startIdx < 0" << "startRow: " << startRow << " numElements: " << numElements << std::endl;
			startIdx = 0;
			numElements = 1;
		}
		MPI_Bcast(&asVec[startIdx], numElements, MPI_DOUBLE, rank, MPI_COMM_WORLD);
		if (rank == 0 && k == n - 1)
		{
			//std::cout << "HELLO: a\n" << a << std::endl;
			a = SparseMatrix<T>(asVec);
		}
		
	}
	
}

template <class T>
void ILU(SparseMatrix<T> & a)
{
	int n = a.getRowCount();
	for (int k = 0; k < n; k++)
	{

		//if (startIdx > 0 && numElements > 0)
		for (int i = k + 1; i < n; i++)
		{
			typename std::multimap<int, std::pair<int, T> >::iterator itKK = a.getIter(k, k);
			//if (itKK != a.data.end())
			//{
			typename std::multimap<int, std::pair<int, T> >::iterator itIK = a.getIter(i, k);
			if (itIK != a.data.end())
			{
				itIK->second.second /= itKK->second.second;
				for (int j = k + 1; j < n; j++)
				{
					typename std::multimap<int, std::pair<int, T> >::iterator itKJ = a.getIter(k, j);
					if (itKJ != a.data.end())
					{
						typename std::multimap<int, std::pair<int, T> >::iterator itIJ = a.getIter(i, j);
						if (itIJ == a.data.end())
						{
							// insertValue checks again
							T val = -1 * itIK->second.second * itKJ->second.second;
							a.data.insert(std::pair<int, std::pair<int, T> >(i, std::make_pair(j, val)));
							std::cout << "adding new pair ("<< i << "," << j << ") of " << val << std::endl;
						}
						else 
						{
							itIJ->second.second -= -1 * (itIK->second.second * itKJ->second.second);
							std::cout  << "modifying ("<< i << "," << j << ")  to " << itIJ->second.second << std::endl;
						}
					}
				}
			}
				//}
		}
			
		//}
		
	}

}