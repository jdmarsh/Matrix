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
	constexpr auto A = Matrix{ {6.0, 2.0}, {1.0, 4.0} };
	static_assert(A.determinant() == 22);

	print_matrix(A);
	
	constexpr auto I = Matrix<float, 3, 3>::Identity();
	print_matrix(I);

	constexpr auto Z = Matrix<float, 3, 3>::Zero();
	print_matrix(Z);
	
	std::cout << "size: " << sizeof(A) << "\n";
	std::cout << "det: " << A.determinant() << "\n";
	std::cout << "trace: " << A.trace() << "\n";
	std::cout << "rank: " << A.rank() << "\n";

	constexpr auto MA = -A;
	print_matrix(MA);
	
	constexpr auto AI = A.inverse();
	print_matrix(AI);
	
	constexpr auto B = A * AI;
	print_matrix(B);
	
	constexpr auto AT = A.transpose();
	print_matrix(AT);
	
	constexpr auto REF = A.REF();
	print_matrix(REF);
	
	constexpr auto RREF = A.RREF();
	print_matrix(RREF);
	
	print_matrix(A);

	print_matrix(A * 2);
	print_matrix(A + A);
	print_matrix(A - A);
	
	auto AC = A;
	AC += AC;
	print_matrix(AC);
	
	AC -= RREF;
	print_matrix(AC);
	
	AC *= 5;
	print_matrix(AC);
}
