#pragma once

#include <array>

template <uint32_t M, uint32_t N>
class Matrix
{
public:
	static Matrix<M, N> Identity()
	requires(M == N)
	{
		Matrix<M, N> identity;
		for (uint32_t i = 0; i < M; ++i)
		{
			for (uint32_t j = 0; j < N; ++j)
			{
				identity[i][j] = (i == j) ? 1 : 0;
			}
		}

		return identity;
	}

	static Matrix<M, N> Zero()
	{
		Matrix<M, N> identity;
		for (uint32_t i = 0; i < M; ++i)
		{
			for (uint32_t j = 0; j < N; ++j)
			{
				identity[i][j] = 0;
			}
		}

		return identity;
	}

	std::array<float, M>::iterator begin()
	{
		return m_matrix.begin();
	}

	std::array<float, M>::const_iterator begin() const
	{
		return m_matrix.begin();
	}

	std::array<float, M>::iterator end()
	{
		return m_matrix.end();
	}

	std::array<float, M>::const_iterator end() const
	{
		return m_matrix.end();
	}

	constexpr float determinant() const
	requires(M > 0 && M == N)
	{
		//Base case
		if constexpr (M == 1)
		{
			return (m_matrix[0][0]);
		}
		else
		{
			//Recursively calculate the determinant
			float determinant = 0;
			for (uint32_t columnIndex = 0; columnIndex < N; ++columnIndex)
			{
				float a = m_matrix[0][columnIndex] * submatrix(0, columnIndex).determinant();
				determinant += ((columnIndex % 2) ? -a : a);
			}

			return determinant;
		}
	}

	float trace() const
	requires(M == N)
	{
		float trace = 0;
		for (unsigned index = 0; index < M; ++index)
		{
			trace += m_matrix[index][index];
		}

		return trace;
	}

	Matrix<M, N> inverse() const
	requires(M > 0 && M == N)
	{
		auto det = determinant();
		if (det == 0)
		{
			return Matrix<M, N>::Zero();
		}

		//If the matrix is a single entry, calculate the inverse differently
		if constexpr (M == 1)
		{
			auto X = Matrix<1, 1>::Identity();
			X *= 1.0f / det;
			return X;
		}

		//Calculate the inverse for larger matrix sizes
		Matrix<M, N> result;
		for (unsigned rowIndex = 0; rowIndex < M; ++rowIndex)
		{
			for (unsigned columnIndex = 0; columnIndex < N; ++columnIndex)
			{
				float a = submatrix(rowIndex, columnIndex).determinant();
				result.m_matrix[rowIndex][columnIndex] = (((rowIndex + columnIndex) % 2) ? -a : a);
			}
		}

		result = result.transpose();
		result *= 1.0f / det;

		return result;
	}

	Matrix<N, M> transpose() const
	requires(M == N)
	{
		Matrix<N, M> result;
		for (unsigned rowIndex = 0; rowIndex < M; ++rowIndex)
		{
			for (unsigned columnIndex = 0; columnIndex < N; ++columnIndex)
			{
				result.m_matrix[columnIndex][rowIndex] = m_matrix[rowIndex][columnIndex];
			}
		}

		return result;
	}

	Matrix<M, N> REF() const
	{
		Matrix<M, N> result;
		result.m_matrix = m_matrix;

		uint32_t x = 0;
		for (uint32_t j = 0; j < N; ++j)
		{
			for (uint32_t i = x; i < M; ++i)
			{
				if (result.m_matrix[i][j])
				{
					result.rowswap(i, x);
					result.rowdivide(x, result.m_matrix[x][j]);
					for (uint32_t a = x + 1; a < M; ++a)
					{
						result.rowminus(a, x, result.m_matrix[a][j]);
					}
					++x;
					break;
				}
			}
		}
		return result;
	}

	Matrix<M, N> RREF() const
	{
		Matrix result = REF();
		for (uint32_t rowIndex = 0; rowIndex < M; ++rowIndex)
		{
			for (uint32_t columnIndex = 0; columnIndex < N; ++columnIndex)
			{
				if (result.m_matrix[rowIndex][columnIndex])
				{
					for (uint32_t a = 0; a < rowIndex; ++a)
					{
						result.rowminus(a, rowIndex, result.m_matrix[a][columnIndex]);
					}
					break;
				}
			}
		}
		return result;
	}

	uint32_t rank() const
	{
		Matrix rref = RREF();
		uint32_t x = 0;
		for (uint32_t rowIndex = 0; rowIndex < M; ++rowIndex)
		{
			for (uint32_t columnIndex = 0; columnIndex < N; ++columnIndex)
			{
				if (rref.m_matrix[rowIndex][columnIndex])
				{
					++x;
					break;
				}
			}
		}
		return x;
	}

	std::array<float, N>& operator[](uint32_t i)
	{
		return m_matrix[i];
	}

	const std::array<float, N>& operator[](uint32_t i) const
	{
		return m_matrix[i];
	}

	void operator+=(const Matrix<M, N>& other)
	{
		for (uint32_t rowIndex = 0; rowIndex < M; ++rowIndex)
		{
			for (uint32_t columnIndex = 0; columnIndex < N; ++columnIndex)
			{
				m_matrix[rowIndex][columnIndex] += other[rowIndex][columnIndex];
			}
		}
	}

	void operator-=(const Matrix<M, N>& other)
	{
		for (uint32_t rowIndex = 0; rowIndex < M; ++rowIndex)
		{
			for (uint32_t columnIndex = 0; columnIndex < N; ++columnIndex)
			{
				m_matrix[rowIndex][columnIndex] -= other[rowIndex][columnIndex];
			}
		}
	}

	template <typename T>
	void operator*=(const T& factor)
	requires (std::is_arithmetic_v<T>)
	{
		for (auto& row : m_matrix)
		{
			for (auto& entry : row)
			{
				entry *= factor;
			}
		}
	}

private:
	std::array<std::array<float, N>, M> m_matrix;

	void rowdivide(uint32_t row, float x)
	{
		for (auto& entry : m_matrix[row])
		{
			entry /= x;
		}
	}

	void rowswap(uint32_t i, uint32_t j)
	{
		m_matrix[i].swap(m_matrix[j]);
	}

	void rowminus(uint32_t i, uint32_t j, float x)
	{
		for (uint32_t a = 0; a < N; ++a)
		{
			m_matrix[i][a] -= (m_matrix[j][a] * x);
		}
	}

	constexpr Matrix<M - 1, N - 1> submatrix(uint32_t i, uint32_t j) const
	{
		Matrix<M - 1, N - 1> result;
		for (uint32_t row = 0; row < M; ++row)
		{
			if (row == i)
			{
				continue;
			}

			auto mappedrow = row > i ? row - 1 : row;
			for (uint32_t col = 0; col < N; ++col)
			{
				if (col == j)
				{
					continue;
				}

				auto mappedcol = col > j ? col - 1 : col;

				result[mappedrow][mappedcol] = m_matrix[row][col];
			}
		}

		return result;
	}
};

template<uint32_t M, uint32_t N, uint32_t P>
Matrix<M, P> operator*(const Matrix<M, N>& mrxA, const Matrix<N, P>& mrxB)
{
	//Create empty matrix of zeros
	Matrix<M, P> result;
	//ToDo: implement begin/end functions to allow users to use range based for loops
	for (uint32_t rowIndex = 0; rowIndex < M; ++rowIndex)
	{
		for (uint32_t columnIndex = 0; columnIndex < P; ++columnIndex)
		{
			result[rowIndex][columnIndex] = 0;
		}
	}

	//Calculate the matrix entries, for loops are rearranged to ensure contiguous data loading
	for (uint32_t rowIndex = 0; rowIndex < M; ++rowIndex)
	{
		for (uint32_t elementIndex = 0; elementIndex < N; ++elementIndex)
		{
			for (uint32_t columnIndex = 0; columnIndex < P; ++columnIndex)
			{
				result[rowIndex][columnIndex] += (mrxA[rowIndex][elementIndex] * mrxB[elementIndex][columnIndex]);
			}
		}
	}

	return result;
}

template<uint32_t M, uint32_t N, typename T>
Matrix<M, N> operator*(const Matrix<M, N>& matrix, const T& factor)
requires (std::is_arithmetic_v<T>)
{
	auto result = matrix;
	for (auto& row : result)
	{
		for (auto& entry : row)
		{
			entry *= factor;
		}
	}
	return result;
}

template<uint32_t M, uint32_t N, typename T>
Matrix<M, N> operator*(const T& factor, const Matrix<M, N>& matrix)
requires (std::is_arithmetic_v<T>)
{
	return matrix * factor;
}

template<uint32_t M, uint32_t N>
Matrix<M, N> operator+(const Matrix<M, N>& mrxA, const Matrix<M, N>& mrxB)
{
	Matrix<M, N> result();
	for (unsigned rowIndex = 0; rowIndex < M; ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < N; ++columnIndex)
		{
			result[rowIndex][columnIndex] = mrxA[rowIndex][columnIndex] + mrxB[rowIndex][columnIndex];
		}
	}

	return result;
}

template<uint32_t M, uint32_t N>
Matrix<M, N> operator-(const Matrix<M, N>& mrxA, const Matrix<M, N>& mrxB)
{
	Matrix<M, N> result();
	for (unsigned rowIndex = 0; rowIndex < M; ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < N; ++columnIndex)
		{
			result[rowIndex][columnIndex] = mrxA[rowIndex][columnIndex] - mrxB[rowIndex][columnIndex];
		}
	}

	return result;
}

template<uint32_t M, uint32_t N>
Matrix<M, N> operator-(const Matrix<M, N>& mrxA)
{
	Matrix<M, N> result();
	for (unsigned rowIndex = 0; rowIndex < M; ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < N; ++columnIndex)
		{
			result[rowIndex][columnIndex] = -mrxA[rowIndex][columnIndex];
		}
	}

	return result;
}

template<uint32_t M, uint32_t N>
bool operator==(const Matrix<M, N>& mrxA, const Matrix<M, N>& mrxB)
{
	for (unsigned rowIndex = 0; rowIndex < M; ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < N; ++columnIndex)
		{
			if (mrxA[rowIndex][columnIndex] != mrxB[rowIndex][columnIndex])
			{
				return false;
			}
		}
	}

	return true;
}

template<uint32_t M, uint32_t N>
bool operator!=(const Matrix<M, N>& mrxA, const Matrix<M, N>& mrxB)
{
	return !(mrxA == mrxB);
}
