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
	const int N = 100;
	const int MAX = 3;
	const int MIN = 0;
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	matrix2D<double> matrixL;
	matrix2D<double> matrixL2;
	SparseMatrix<double> sm;
	double startP, endP, startS, endS, startSadd, endSadd, startPadd, endPadd;
	if (rank == 0)
	{
		matrixL = matrix2D<double>(N, N);
		//matrixL.makePositiveDefinite();
		//matrixL.fillRandomlyLowerTriangular(MIN, MAX);
		matrixL2 = matrix2D<double>(N, N, 1.0);
		// std::cout << "matrixL: \n" << matrixL << std::endl;
		// std::cout << "matrixL2: \n" << matrixL2 << std::endl;
		startS = MPI_Wtime();
		//multiply(matrixL, matrixL2);
		endS = MPI_Wtime();
		startSadd = MPI_Wtime();
		add(matrixL, matrixL2);
		endSadd = MPI_Wtime();

		//std::cout << "submatrixL2: \n" << matrixL.createSubMatrix(0, 3, 0, 1) << std::endl;
		//matrix2D<double> product = multiply(matrixL, matrixL2);
		//std::cout << "multiply done: " << std::endl;
		//std::cout << "product: \n" << product << std::endl;

	}

	startP = MPI_Wtime();
	parallelMult(rank, size, matrixL, matrixL2);
	endP = MPI_Wtime();

	startPadd = MPI_Wtime();
	parallelAdd(rank, size, matrixL, matrixL2);
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
