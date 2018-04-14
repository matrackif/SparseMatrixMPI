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
	matrix2D<double> matrixL;
	matrix2D<double> matrixL2;
	matrix2D<double> mtx;
	//SparseMatrix<double> sm;
	/*
	sm.insertValue(3, 2, 6);
	sm.insertValue(3, 1, 3);
	sm.insertValue(3, 2, 4);
	sm.insertValue(0, 4, 5.3);
	sm.insertValue(1, 2, 11);
	sm.insertValue(1, 2, 11);
	sm.insertValue(1, 2, 0.0);
	*/
	
	if (rank == 0)
	{
		SparseMatrix<double> sm = SparseMatrix<double>("../../test.txt");
		SparseMatrix<double> sm2 = SparseMatrix<double>("../../test2.txt");
		std::cout << "SM: \n" << sm << std::endl;
		std::cout << "RowCount: " << sm.getRowCount() << " ColumnCount: " << sm.getColumnCount() 
			<< " NonZeroCount: " << sm.getNumNonZeroElements() << std::endl;

		std::cout << "SM2: \n" << sm2 << std::endl;
		std::cout << "RowCount: " << sm2.getRowCount() << " ColumnCount: " << sm2.getColumnCount()
			<< " NonZeroCount: " << sm2.getNumNonZeroElements() << std::endl;

		SparseMatrix<double> res = sm.multiply(sm2);
		std::cout << "res: \n" << res << std::endl;
		std::cout << "RowCount: " << res.getRowCount() << " ColumnCount: " << res.getColumnCount()
			<< " NonZeroCount: " << res.getNumNonZeroElements() << std::endl;

		res = sm.add(res);
		std::cout << "res after adding: \n" << res << std::endl;
		std::cout << "RowCount: " << res.getRowCount() << " ColumnCount: " << res.getColumnCount()
			<< " NonZeroCount: " << res.getNumNonZeroElements() << std::endl;
	}
	
	
	double startP, endP, startS, endS, startSadd, endSadd, startPadd, endPadd;
	if (rank == 0)
	{
		// matrixL = matrix2D<double>(N, N);
		//matrixL.makePositiveDefinite();
		//matrixL.fillRandomlyLowerTriangular(MIN, MAX);
		// matrixL2 = matrix2D<double>(N, N, 1.0);
		// std::cout << "matrixL: \n" << matrixL << std::endl;
		// std::cout << "matrixL2: \n" << matrixL2 << std::endl;
		//mtx = matrix2D<double>("../../cavity08.mtx");
		//std::cout << "mtx: \n" << mtx << std::endl;
		//std::cout << "RowCount: " << mtx.getRowCount() << " ColumnCount: " << mtx.getColumnCount() << std::endl;
		startS = MPI_Wtime();
		//multiply(matrixL, matrixL2);
		//multiply(mtx, mtx);
		endS = MPI_Wtime();
		startSadd = MPI_Wtime();
		//add(matrixL, matrixL2);
		//add(mtx, mtx);
		endSadd = MPI_Wtime();

		//std::cout << "submatrixL2: \n" << matrixL.createSubMatrix(0, 3, 0, 1) << std::endl;
		//matrix2D<double> product = multiply(matrixL, matrixL2);
		//std::cout << "multiply done: " << std::endl;
		//std::cout << "product: \n" << product << std::endl;

	}

	startP = MPI_Wtime();
	//parallelMult(rank, size, mtx, mtx);
	endP = MPI_Wtime();

	startPadd = MPI_Wtime();
	//parallelAdd(rank, size, mtx, mtx);
	endPadd = MPI_Wtime();

	if (rank == 0)
	{
		std::cout << "Sequential mult time: " << endS - startS << "s" << std::endl;
		std::cout << "Parallel mult time: " << endP - startP << "s" << std::endl;
		std::cout << "Sequential add time: " << endSadd - startSadd << "s" << std::endl;
		std::cout << "Parallel add time: " << endPadd - startPadd << "s" << std::endl;
	}
	

	MPI_Finalize();
}
