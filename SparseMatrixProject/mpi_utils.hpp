#pragma once
#include "matrix.hpp"
#include "mpi.h"


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