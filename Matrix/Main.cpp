#include <iostream>
#include "Matrix.h"

int main()
{
	Matrix A(2, 2);
	A[0][0] = 6;
	A[0][1] = 2;
	A[1][0] = 0;
	A[1][1] = 4;
	try
	{
		Matrix E = A * A;
		std::cout << E[0][0] << ",";
		std::cout << E[0][1] << std::endl;
		std::cout << E[1][0] << ",";
		std::cout << E[1][1];
	}
	catch (std::string error)
	{
		std::cout << error << std::endl;
	}
	return 0;
}
