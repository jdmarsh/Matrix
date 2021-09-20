#include "Matrix.h"

Matrix::Matrix(unsigned i, unsigned j)
{
	matrix.resize(i);
	for (std::vector<float>& row : matrix)
	{
		row.resize(j);
	}
}

Matrix::Matrix(std::initializer_list<std::initializer_list<float>> list)
{
	//Set up the vectors for the matrix
	unsigned rowCount = list.size();
	if (rowCount > 0)
	{
		matrix.resize(rowCount);
		unsigned index = 0;
		for (std::initializer_list<float> row : list)
		{
			matrix[index].resize(row.size());
			++index;
		}
	}
	//Fill matrix with data
	unsigned rowIndex = 0;
	for (std::initializer_list<float> row : list)
	{
		unsigned columnIndex = 0;
		for (float value : row)
		{
			matrix[rowIndex][columnIndex] = value;
			++columnIndex;
		}
		++rowIndex;
	}
}

unsigned Matrix::cols()
{
	if (matrix.size() == 0)
	{
		return 0;
	}
	return matrix[0].size();
}

unsigned Matrix::rows()
{
	return matrix.size();
}

float Matrix::determinant()
{
	if (rows() == 0 || cols() == 0)
	{
		throw std::string("Unable to calculate the deteminant of an empty matrix");
	}
	if (rows() != cols())
	{
		throw std::string("Unable to calculate the determinant of a non square matrix");
	}

	//Base case
	if (rows() == 1)
	{
		return (matrix[0][0]);
	}

	//Recursively calculate the determinant
	float determinant = 0;
	for (unsigned columnIndex = 0; columnIndex < cols(); ++columnIndex)
	{
		float a = matrix[0][columnIndex] * newmrx(0, columnIndex).determinant();
		determinant += ((columnIndex % 2) ? -a : a);
	}

	return determinant;
}

float Matrix::trace()
{
	if (rows() != cols())
	{
		throw std::string("Unable to calculate the trace of a non square matrix");
	}

	float trace = 0;
	for (unsigned index = 0; index < rows(); ++index)
	{
		trace += matrix[index][index];
	}

	return trace;
}

Matrix Matrix::inverse()
{
	if (rows() == 0 || cols() == 0)
	{
		throw std::string("Unable to calculate the inverse of an empty matrix");
	}
	if (rows() != cols())
	{
		throw std::string("Unable to calculate the inverse of a non square matrix");
	}

	float det = determinant();
	if (det == 0)
	{
		throw std::string("Unable to calculate the inverse of a matrix with a determinant of 0");
	}

	//If the matrix is a single entry, calculate the inverse differently
	if (matrix.size() == 1)
	{
		Matrix X = id(1);
		X.scale(1.0f / det);
		return X;
	}

	//Calculate the inverse for larger matrix sizes
	Matrix A(rows(), cols());
	for (unsigned rowIndex = 0; rowIndex < rows(); ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < cols(); ++columnIndex)
		{
			float a = newmrx(rowIndex, columnIndex).determinant();
			A.matrix[rowIndex][columnIndex] = (((rowIndex + columnIndex) % 2) ? -a : a);
		}
	}

	A = A.transpose();
	A.scale(1.0f / det);

	return A;
}

Matrix Matrix::transpose()
{
	if (rows() == 0 || cols() == 0)
	{
		throw std::string("Unable to calculate the inverse of an empty matrix");
	}

	Matrix B(cols(), rows());
	for (unsigned rowIndex = 0; rowIndex < rows(); ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < cols(); ++columnIndex)
		{
			B.matrix[columnIndex][rowIndex] = matrix[rowIndex][columnIndex];
		}
	}

	return B;
}

Matrix Matrix::REF()
{
	Matrix M(rows(), cols());
	M.matrix = matrix;
	unsigned x = 0;
	for (unsigned j = 0; j < cols(); ++j)
	{
		for (unsigned i = x; i < rows(); ++i)
		{
			if (M.matrix[i][j])
			{
				M.rowswap(i, x);
				M.rowdivide(x, M.matrix[x][j]);
				for (unsigned a = x + 1; a < rows(); ++a)
				{
					M.rowminus(a, x, M.matrix[a][j]);
				}
				x++;
				break;
			}
		}
	}
	return M;
}

Matrix Matrix::RREF()
{
	Matrix M = REF();
	for (unsigned rowIndex = 0; rowIndex < M.rows(); ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < M.cols(); ++columnIndex)
		{
			if (M.matrix[rowIndex][columnIndex])
			{
				for (unsigned a = 0; a < rowIndex; ++a)
				{
					M.rowminus(a, rowIndex, M.matrix[a][columnIndex]);
				}
				break;
			}
		}
	}
	return M;
}

unsigned Matrix::rank()
{
	Matrix M = RREF();
	unsigned x = 0;
	for (unsigned rowIndex = 0; rowIndex < M.rows(); ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < M.cols(); ++columnIndex)
		{
			if (M.matrix[rowIndex][columnIndex])
			{
				x++;
				break;
			}
		}
	}
	return x;
}

void Matrix::scale(float x)
{
	for (unsigned rowIndex = 0; rowIndex < rows(); ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < cols(); ++columnIndex)
		{
			matrix[rowIndex][columnIndex] *= x;
		}
	}
}

std::vector<float>& Matrix::operator[](unsigned i)
{
	return matrix[i];
}

void Matrix::operator+=(Matrix& A)
{
	if ((rows() != A.rows()) || (cols() != A.cols()))
	{
		throw std::string("Matrix sizes do not align correctly for addition");
	}

	for (unsigned rowIndex = 0; rowIndex < rows(); ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < cols(); ++columnIndex)
		{
			matrix[rowIndex][columnIndex] += A[rowIndex][columnIndex];
		}
	}
}

void Matrix::operator-=(Matrix& A)
{
	if ((rows() != A.rows()) || (cols() != A.cols()))
	{
		throw std::string("Matrix sizes do not align correctly for subtraction");
	}

	for (unsigned rowIndex = 0; rowIndex < rows(); ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < cols(); ++columnIndex)
		{
			matrix[rowIndex][columnIndex] -= A[rowIndex][columnIndex];
		}
	}
}

void Matrix::rowdivide(unsigned i, float x)
{
	for (unsigned a = 0; a < cols(); ++a)
	{
		matrix[i][a] /= x;
	}
}

void Matrix::rowswap(unsigned i, unsigned j)
{
	std::vector<float> vec;
	vec = matrix[i];
	matrix[i] = matrix[j];
	matrix[j] = vec;
}

void Matrix::rowminus(unsigned i, unsigned j, float x)
{
	for (unsigned a = 0; a < cols(); ++a)
	{
		matrix[i][a] -= (matrix[j][a] * x);
	}
}

Matrix Matrix::newmrx(unsigned i, unsigned j)
{
	if (rows() == 0) return id(0);
	Matrix mrx(rows(), cols());
	mrx.matrix = matrix;
	for (unsigned a = 0; a < rows(); ++a)
	{
		mrx.matrix[a].erase(mrx.matrix[a].begin() + j);
	}
	mrx.matrix.erase(mrx.matrix.begin() + i);
	return mrx;
}

Matrix id(unsigned dimention)
{
	Matrix identity(dimention, dimention);
	for (unsigned index = 0; index < dimention; ++index)
	{
		identity[index][index] = 1;
	}

	return identity;
}

Matrix operator*(Matrix& A, Matrix& B)
{
	if (A.cols() != B.rows())
	{
		throw std::string("Matrix sizes do not align correctly for multiplication");
	}

	//Create empty matrix of zeros
	Matrix result(A.rows(), B.cols());
	//ToDo: implement begin/end functions to allow users to use range based for loops
	for (unsigned rowIndex = 0; rowIndex < result.rows(); ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < result.cols(); ++columnIndex)
		{
			result[rowIndex][columnIndex] = 0;
		}
	}

	//Calculate the matrix entries, for loops are rearranged to ensure contiguous data loading
	for (unsigned rowIndex = 0; rowIndex < result.rows(); ++rowIndex)
	{
		for (unsigned elementIndex = 0; elementIndex < A.cols(); ++elementIndex)
		{
			for (unsigned columnIndex = 0; columnIndex < result.cols(); ++columnIndex)
			{
				result[rowIndex][columnIndex] += (A[rowIndex][elementIndex] * B[elementIndex][columnIndex]);
			}
		}
	}

	return result;
}

Matrix operator+(Matrix& A, Matrix& B)
{
	if ((A.rows() != B.rows()) || (A.cols() != B.cols()))
	{
		throw std::string("Matrix sizes do not align correctly for addition");
	}

	Matrix result(A.rows(), A.cols());
	for (unsigned rowIndex = 0; rowIndex < A.rows(); ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < A.cols(); ++columnIndex)
		{
			result[rowIndex][columnIndex] = A[rowIndex][columnIndex] + B[rowIndex][columnIndex];
		}
	}

	return result;
}

Matrix operator-(Matrix& A, Matrix& B)
{
	if ((A.rows() != B.rows()) || (A.cols() != B.cols()))
	{
		throw std::string("Matrix sizes do not align correctly for subtraction");
	}

	Matrix result(A.rows(), A.cols());
	for (unsigned rowIndex = 0; rowIndex < A.rows(); ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < A.cols(); ++columnIndex)
		{
			result[rowIndex][columnIndex] = A[rowIndex][columnIndex] - B[rowIndex][columnIndex];
		}
	}

	return result;
}

Matrix operator-(Matrix& A)
{
	Matrix result(A.rows(), A.cols());
	for (unsigned rowIndex = 0; rowIndex < A.rows(); ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < A.cols(); ++columnIndex)
		{
			result[rowIndex][columnIndex] = -A[rowIndex][columnIndex];
		}
	}

	return result;
}

bool operator==(Matrix& A, Matrix& B)
{
	if ((A.rows() != B.rows()) || (A.cols() != B.cols()))
	{
		return false;
	}

	for (unsigned rowIndex = 0; rowIndex < A.rows(); ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < A.cols(); ++columnIndex)
		{
			if (A[rowIndex][columnIndex] != B[rowIndex][columnIndex])
			{
				return false;
			}
		}
	}

	return true;
}

bool operator!=(Matrix& A, Matrix& B)
{
	return !(A == B);
}
