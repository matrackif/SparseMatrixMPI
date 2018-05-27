#include <mpi.h>
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <cmath>
#include "matrix.hpp"
#include "sparse_matrix.hpp"
#include "mpi_utils.hpp"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	const int N = 5;
	const int MAX = 3;
	const int MIN = 0;
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	SparseMatrix<double> testA, testX, testB, matrixFromFile;
	SparseMatrix<double> lowerTriangSparse, lowerTriangSparseTransposed, testPosDef;
	if (rank == 0)
	{
		lowerTriangSparse = SparseMatrix<double>(7, 7);
		lowerTriangSparse.fillRandomlyLowerTriangular(1, 3, 1);
		lowerTriangSparseTransposed = transpose(lowerTriangSparse);
		testPosDef = lowerTriangSparse.multiply(lowerTriangSparseTransposed);
		std::cout << "lowerTriangSparse: \n" << lowerTriangSparse << std::endl;
		std::cout << "lowerTriangSparseTransposed: \n" << lowerTriangSparseTransposed << std::endl;
		std::cout << "testPosDef: \n" << testPosDef << std::endl;
		testA = SparseMatrix<double>(3, 3, 0);
		testX = SparseMatrix<double>(3, 1, 0);
		testB = SparseMatrix<double>(3, 1, 0);
		testA.insertValue(0, 0, 1);
		testA.insertValue(1, 0, 1);
		testA.insertValue(1, 1, 1);
		testA.insertValue(2, 0, 1);
		testA.insertValue(2, 1, 1);
		testA.insertValue(2, 2, 1);
		SparseMatrix<double> testATransposed = transpose(testA);
		testA = testA.multiply(testATransposed);
		testB.insertValue(0, 0, 3);
		testB.insertValue(1, 0, 5);
		testB.insertValue(2, 0, 6);
		std::cout << "testA: \n" << testA << std::endl;
		std::cout << "RowCount: " << testA.getRowCount() << " ColumnCount: " << testA.getColumnCount()
			<< " NonZeroCount: " << testA.getNumNonZeroElements() << std::endl;
		std::cout << "testX: \n" << testX << std::endl;
		std::cout << "RowCount: " << testX.getRowCount() << " ColumnCount: " << testX.getColumnCount()
			<< " NonZeroCount: " << testX.getNumNonZeroElements() << std::endl;
		std::cout << "testB: \n" << testB << std::endl;
		std::cout << "RowCount: " << testB.getRowCount() << " ColumnCount: " << testB.getColumnCount()
			<< " NonZeroCount: " << testB.getNumNonZeroElements() << std::endl;
		/*
		matrixFromFile = SparseMatrix<double>("../../bcsstk05.mtx");
		std::cout << "Matrix from file before decomp:\n" << matrixFromFile << "\n";
		incompleteCholeskyDecomp(matrixFromFile);
		std::cout << "Matrix from file after decomp:\n" << matrixFromFile << "\n";
		*/
	}
	conjugateGradient(rank, size, testA, testB, testX);
	
	//incompleteCholeskyDecomp(testPosDef);
	parallelILU(rank, size, testPosDef);
	
	if (rank == 0)
	{
		/*
		std::cout << "testPosDef was factorized into matrix L:\n" << testPosDef << std::endl;
		SparseMatrix<double> factorizedPosDefTransposed = transpose(testPosDef);
		SparseMatrix<double> multResult = testPosDef.multiply(factorizedPosDefTransposed);
		std::cout << "L * LT:\n" << multResult << std::endl;
		*/
		std::cout << "Solution X found for matrix testA:\n" << testX << std::endl;
		std::cout << "parallelILU() decomposition found:\n" << testPosDef << std::endl;
	}
	
	MPI_Finalize();
}
