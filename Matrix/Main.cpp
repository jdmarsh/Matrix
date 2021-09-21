#include <iostream>
#include "Matrix.h"

template<std::floating_point T, uint32_t M, uint32_t N>
void print_matrix(const Matrix<T, M, N>& matrix)
{
	for (uint32_t rowIndex = 0; rowIndex < M; ++rowIndex)
	{
		if (rowIndex > 0)
		{
			std::cout << "\n";
		}

		for (uint32_t columnIndex = 0; columnIndex < N; ++columnIndex)
		{
			if (columnIndex > 0)
			{
				std::cout << ",";
			}
			std::cout << matrix[rowIndex][columnIndex];
		}
	}
	std::cout << "\n";
}

int main()
{
	auto A = Matrix{ {6.0, 2.0}, {0.0, 4.0} };
	
	print_matrix(A);
	
	constexpr auto I = Matrix<float, 3, 3>::Identity();
	
	print_matrix(I);
	auto Z = Matrix<float, 3, 3>::Zero();
	print_matrix(Z);
	
	std::cout << "det: " << A.determinant() << "\n";
	std::cout << "trace: " << A.trace() << "\n";
	std::cout << "rank: " << A.rank() << "\n";
	
	auto AI = A.inverse();
	print_matrix(AI);
	
	auto B = A * AI;
	print_matrix(B);
	
	auto AT = A.transpose();
	print_matrix(AT);
	
	auto REF = A.REF();
	print_matrix(REF);
	
	auto RREF = A.RREF();
	print_matrix(RREF);
	
	print_matrix(A);
	
	A += A;
	print_matrix(A);
	
	A -= RREF;
	print_matrix(A);
	
	A *= 5;
	print_matrix(A);
}
