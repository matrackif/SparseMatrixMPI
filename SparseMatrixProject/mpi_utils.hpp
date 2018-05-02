#pragma once
#include "matrix.hpp"
#include "mpi.h"

template <typename T>
void concat(std::vector<T>& a, const std::vector<T>& b)
{
	a.reserve(a.size() + b.size());
	a.insert(a.end(), b.begin(), b.end());
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
void parallelAdd(int rank, int size, SparseMatrix<T> & m1, SparseMatrix<T> & m2)
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
		return;
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
		std::cout << "parallelAdd() Rank=" << rank << " numRowsPerProcess: " << numRowsPerProcess << std::endl;
		std::cout << "parallelAdd() Rank=" << rank << " numRowsPerProcessM2: " << numRowsPerProcessM2 << std::endl;
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
		std::cout << "SparseMatrix parallelAdd() Finished! m1: \n" << m1 << "m2: \n" << m2 << " m1 + m2: \n" << result << std::endl;
	}
	else
	{
		std::cout << "parallelAdd() Rank=" << rank << " numRowsPerProcess: " << numRowsPerProcess << std::endl;
		std::cout << "parallelAdd() Rank=" << rank << " numRowsPerProcessM2: " << numRowsPerProcessM2 << std::endl;
		int m1Size, m2Size;
		MPI_Recv(&m1Size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&m2Size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		std::cout << "parallelAdd() Rank=" << rank << " m1Size: " << m1Size << " m2Size: " << m2Size << std::endl;
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
}

template <class T>
void parallelMult(int rank, int size, SparseMatrix<T> & m1, SparseMatrix<T> & m2)
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
		std::cout << "parallelMult() Rank=" << rank << " numRowsPerProcess: " << numRowsPerProcess << std::endl;
		//std::cout << "matrix m1: \n" << m1 << std::endl;
		//std::cout << "matrix m2: \n" << m2 << std::endl;
		result = SparseMatrix<T>(m1.getRowCount(), m2.getColumnCount());
		std::vector<T> resultAsVec = result.toVector();
		//std::cout << "resultAsVec: \n";
		for (int k = 0; k < resultAsVec.size(); k++)
		{
			std::cout << resultAsVec[k] << ",";
		}
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
		std::cout << "SparseMatrix parallelMult() Finished! m1: \n" << m1 << "m2: \n" << m2 << " m1 * m2: \n" << result << std::endl;
	}
	else
	{
		std::cout << "parallelMult() Rank=" << rank << " numRowsPerProcess: " << numRowsPerProcess << std::endl;

		int m1Size, m2Size;
		MPI_Recv(&m1Size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&m2Size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		std::cout << "parallelMult() Rank=" << rank << " m1Size: " << m1Size << " m2Size: " << m2Size << std::endl;
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
}

template <class T>
void conjugateGradient(SparseMatrix<T> A, SparseMatrix<T> b, SparseMatrix<T> x)
{
	SparseMatrix<T> r = 
}
