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
	SparseMatrix<double> sm, sm2, sm3, smk;
    matrix2D<double> m;
	
	double startP, endP, startS, endS, startSadd, endSadd, startPadd, endPadd;
    double startPmultSparse, endPmultSparse, startPaddSparse, endPaddSparse, startSaddSparse, endSaddSparse, startSmultSparse, endSmultSparse;
	if (rank == 0)
	{
		
		sm = SparseMatrix<double>("../../test1.txt");
		sm2 = SparseMatrix<double>("../../test2.txt");
		sm3 = SparseMatrix<double>("../../bcsstk05.mtx");
		smk = SparseMatrix<double>(10, 10);
		smk.makeKdiagonal(5);
		m = matrix2D<double>("../../bcsstk05.mtx");
		std::cout << "SMK: \n" << smk << std::endl;
		std::cout << "RowCount: " << smk.getRowCount() << " ColumnCount: " << smk.getColumnCount()
			<< " NonZeroCount: " << smk.getNumNonZeroElements() << std::endl;

		std::vector<std::vector<double> > vec = sm2.getUniqueRows();
		for (int i = 0; i < vec.size(); i++)
		{
			for (int j = 0; j < vec[i].size(); j++)
			{
				std::cout << vec[i][j] << ",";
			}
			std::cout << std::endl;
		}

		std::cout << "SM: \n" << sm << std::endl;
		std::cout << "RowCount: " << sm.getRowCount() << " ColumnCount: " << sm.getColumnCount()
			<< " NonZeroCount: " << sm.getNumNonZeroElements() << std::endl;

		std::cout << "SM2: \n" << sm2 << std::endl;
		std::cout << "RowCount: " << sm2.getRowCount() << " ColumnCount: " << sm2.getColumnCount()
			<< " NonZeroCount: " << sm2.getNumNonZeroElements() << std::endl;

		SparseMatrix<double> res = sm.multiply(sm2);
		std::cout << "res after multiplying sm * sm2: \n" << res << std::endl;
		std::cout << "RowCount: " << res.getRowCount() << " ColumnCount: " << res.getColumnCount()
			<< " NonZeroCount: " << res.getNumNonZeroElements() << std::endl;

		res = sm.add(sm2);
		std::cout << "res after adding sm + sm2: \n" << res << std::endl;
		std::cout << "RowCount: " << res.getRowCount() << " ColumnCount: " << res.getColumnCount()
			<< " NonZeroCount: " << res.getNumNonZeroElements() << std::endl;

		startS = MPI_Wtime();
		multiply(m, m);
		endS = MPI_Wtime();

		startSadd = MPI_Wtime();
		add(m, m);
		endSadd = MPI_Wtime();

        startSmultSparse = MPI_Wtime();
        sm3.multiply(sm3);
        endSmultSparse = MPI_Wtime();

        startSaddSparse = MPI_Wtime();
        sm3.add(sm3);
        endSaddSparse = MPI_Wtime();

	}

	startP = MPI_Wtime();
	parallelMult(rank, size, m, m);
	endP = MPI_Wtime();

	startPadd = MPI_Wtime();
    parallelAdd(rank, size, m, m);
	endPadd = MPI_Wtime();

    startPmultSparse = MPI_Wtime();
	parallelMult(rank, size, sm3, sm3);
	endPmultSparse = MPI_Wtime();

	startPaddSparse = MPI_Wtime();
    parallelAdd(rank, size, sm3, sm3);
	endPaddSparse = MPI_Wtime();

	if (rank == 0)
	{
		std::cout << "Sequential mult time: " << endS - startS << "s" << std::endl;
		std::cout << "Parallel mult time: " << endP - startP << "s" << std::endl;
		std::cout << "Sequential add time: " << endSadd - startSadd << "s" << std::endl;
		std::cout << "Parallel add time: " << endPadd - startPadd << "s" << std::endl;

        std::cout << "Sparse sequential mult time: " << endSmultSparse - startSmultSparse << "s" << std::endl;
		std::cout << "Sparse parallel mult time: " << endPmultSparse - startPmultSparse << "s" << std::endl;
		std::cout << "Sparse sequential add time: " << endSaddSparse - startSaddSparse << "s" << std::endl;
		std::cout << "Sparse parallel add time: " << endPaddSparse - startPaddSparse << "s" << std::endl;
	}
	

	MPI_Finalize();
}
