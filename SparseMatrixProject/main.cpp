#include <mpi.h>
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <cmath>
#include "matrix.hpp"
#include "mpi_utils.hpp"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	const int N = 6;
	const int MAX = 10;
	const int MIN = 0;
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	matrix2D<double> matrixL;
	matrix2D<double> matrixL2;
	if (rank == 0)
	{
		matrixL = matrix2D<double>(N, N);
		//matrixL.makePositiveDefinite();
		matrixL.fillRandomlyLowerTriangular(0, 3);
		matrixL2 = matrix2D<double>(N, N, 1.0);
		std::cout << "matrixL: \n" << matrixL << std::endl;
		std::cout << "matrixL2: \n" << matrixL2 << std::endl;
		
		std::cout << "submatrixL2: \n" << matrixL2.createSubMatrix(0, 3, 0, 2) << std::endl;
		/*
		matrix2D<double> product = multiply(matrixL, matrixL2);
		std::cout << "multiply done: " << std::endl;
		std::cout << "product: \n" << product << std::endl;
		*/
	}
	parallelMult(rank, size, matrixL, matrixL2);
	parallelAdd(rank, size, matrixL, matrixL2);
	MPI_Finalize();
	//system("pause");
}


/*
if (rank == 0)
{
for (int i = 0; i < N; i++)
{
for (int j = 0; j < N; j++)
{
if (i >= j)
{
double randNum = (double)(MIN + (rand() % static_cast<int>(MAX - MIN + 1)));
randNum = (randNum == 0.0) ? 1.0 : randNum;
matrixL(i, j) = randNum;
}
else
{
matrixL(i, j) = 0.0;
}
}
}
std::cout << "matrixL:" << std::endl << matrixL << std::endl;
matrix2D<double> matrixLT = transpose(matrixL);
std::cout << "matrixL transposed:" << std::endl << matrixLT << std::endl;
matrixL = multiply(matrixL, matrixLT);
}
MPI_Bcast(matrixL.getData(0), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD); // rank == 0
if (rank == 0)
{
std::cout << "matrix to decompose: " << std::endl << matrixL << std::endl;
}
// std::cout << "Process with rank " << rank << " sees matrixL as " << matrixL << std::endl;

for (int j = 0; j < N; j++)
{
if (rank == 0)
{
for (int i = 0; i < j; i++)
{
matrixL(i, j) = 0.0;
}
}
if (j % size == rank)
{
for (int k = 0; k < j; k++)
{
matrixL(j, j) -= matrixL(j, k) * matrixL(j, k);
}
matrixL(j, j) = sqrt(matrixL(j, j));

}
MPI_Bcast(matrixL.getData(j), N, MPI_DOUBLE, j % size, MPI_COMM_WORLD);
for (int i = j + 1; i < N; i++)
{
if (i % size == rank)
{
for (int k = 0; k < j; k++)
{
matrixL(i, j) -= matrixL(i, k) * matrixL(j, k);
}
matrixL(i, j) /= matrixL(j, j);
}
}
}

MPI_Barrier(MPI_COMM_WORLD);
if (rank == 0)
{
std::cout << "Final matrixL" << std::endl << matrixL << std::endl;
matrix2D<double> matrixLT = transpose(matrixL);
std::cout << "Final matrixTransposed" << std::endl << matrixLT << std::endl;
std::cout << "L times Lt" << std::endl << multiply(matrixL, matrixLT) << std::endl;
}
*/