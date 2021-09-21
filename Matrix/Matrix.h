#pragma once

#include <array>

template <std::floating_point T, uint32_t M, uint32_t N>
struct Matrix
{
	static auto constexpr rows = M;
	static auto constexpr cols = N;
	static auto constexpr size = M * N;

private:
	template <typename Input, std::size_t... Is>
	auto constexpr init(Input&& in, std::index_sequence<Is...>)
	{
		return std::array<T, sizeof...(Is)>{ static_cast<T>(in[Is])... };
	}

public:
	constexpr Matrix() = default;

	template <typename... Ts>
	constexpr Matrix(Ts(&&...in)[N])
	requires (sizeof...(Ts) == M && ((std::is_nothrow_convertible_v<Ts, T>) && ...))
	: m_matrix{ init(std::move(in), std::make_index_sequence<N>{})... }
	{
	}

	constexpr static Matrix<T, M, N> Identity()
	requires(M == N)
	{
		Matrix<T, M, N> identity;
		for (uint32_t i = 0; i < M; ++i)
		{
			for (uint32_t j = 0; j < N; ++j)
			{
				identity[i][j] = static_cast<T>(i == j);
			}
		}

		return identity;
	}

	constexpr static Matrix<T, M, N> Zero()
	{
		Matrix<T, M, N> identity;
		for (uint32_t i = 0; i < M; ++i)
		{
			for (uint32_t j = 0; j < N; ++j)
			{
				identity[i][j] = static_cast<T>(0);
			}
		}

		return identity;
	}

	constexpr std::array<T, M>::iterator begin()
	{
		return m_matrix.begin();
	}

	constexpr std::array<T, M>::const_iterator begin() const
	{
		return m_matrix.begin();
	}

	constexpr std::array<T, M>::iterator end()
	{
		return m_matrix.end();
	}

	constexpr std::array<T, M>::const_iterator end() const
	{
		return m_matrix.end();
	}

	constexpr T determinant() const
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
			T determinant = 0;
			for (uint32_t columnIndex = 0; columnIndex < N; ++columnIndex)
			{
				T a = m_matrix[0][columnIndex] * submatrix(0, columnIndex).determinant();
				determinant += ((columnIndex % 2) ? -a : a);
			}

			return determinant;
		}
	}

	constexpr T trace() const
	requires(M == N)
	{
		T trace = 0;
		for (unsigned index = 0; index < M; ++index)
		{
			trace += m_matrix[index][index];
		}

		return trace;
	}

	constexpr Matrix<T, M, N> inverse() const
	requires(M > 0 && M == N)
	{
		auto det = determinant();
		if (det == 0)
		{
			return Matrix<T, M, N>::Zero();
		}

		//If the matrix is a single entry, calculate the inverse differently
		if constexpr (M == 1)
		{
			auto X = Matrix<T, 1, 1>::Identity();
			X *= 1.0f / det;
			return X;
		}

		//Calculate the inverse for larger matrix sizes
		Matrix<T, M, N> result;
		for (unsigned rowIndex = 0; rowIndex < M; ++rowIndex)
		{
			for (unsigned columnIndex = 0; columnIndex < N; ++columnIndex)
			{
				T a = submatrix(rowIndex, columnIndex).determinant();
				result.m_matrix[rowIndex][columnIndex] = (((rowIndex + columnIndex) % 2) ? -a : a);
			}
		}

		result = result.transpose();
		result *= static_cast<T>(1) / det;

		return result;
	}

	constexpr Matrix<T, N, M> transpose() const
	requires(M == N)
	{
		Matrix<T, N, M> result;
		for (unsigned rowIndex = 0; rowIndex < M; ++rowIndex)
		{
			for (unsigned columnIndex = 0; columnIndex < N; ++columnIndex)
			{
				result.m_matrix[columnIndex][rowIndex] = m_matrix[rowIndex][columnIndex];
			}
		}

		return result;
	}

	constexpr Matrix<T, M, N> REF() const
	{
		Matrix<T, M, N> result;
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

	constexpr Matrix<T, M, N> RREF() const
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

	constexpr uint32_t rank() const
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

	constexpr std::array<T, N>& operator[](uint32_t i)
	{
		return m_matrix[i];
	}

	constexpr const std::array<T, N>& operator[](uint32_t i) const
	{
		return m_matrix[i];
	}

	template<std::floating_point TO>
	constexpr void operator+=(const Matrix<TO, M, N>& other)
	{
		for (uint32_t rowIndex = 0; rowIndex < M; ++rowIndex)
		{
			for (uint32_t columnIndex = 0; columnIndex < N; ++columnIndex)
			{
				m_matrix[rowIndex][columnIndex] += other[rowIndex][columnIndex];
			}
		}
	}

	template<std::floating_point TO>
	constexpr void operator-=(const Matrix<TO, M, N>& other)
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
	constexpr void operator*=(const T& factor)
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
	std::array<std::array<T, N>, M> m_matrix;

	constexpr void rowdivide(uint32_t row, T x)
	{
		for (auto& entry : m_matrix[row])
		{
			entry /= x;
		}
	}

	constexpr void rowswap(uint32_t i, uint32_t j)
	{
		m_matrix[i].swap(m_matrix[j]);
	}

	constexpr void rowminus(uint32_t i, uint32_t j, T x)
	{
		for (uint32_t a = 0; a < N; ++a)
		{
			m_matrix[i][a] -= (m_matrix[j][a] * x);
		}
	}

	constexpr Matrix<T, M - 1, N - 1> submatrix(uint32_t i, uint32_t j) const
	{
		Matrix<T, M - 1, N - 1> result;
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

template <class... Ts, std::size_t N>
Matrix(Ts(&&..._)[N])->Matrix<std::common_type_t<Ts...>, sizeof...(Ts), N>;

template<std::floating_point TA, std::floating_point TB, uint32_t M, uint32_t N, uint32_t P>
Matrix<std::common_type_t<TA, TB>, M, P> operator*(const Matrix<TA, M, N>& mrxA, const Matrix<TB, N, P>& mrxB)
{
	//Create empty matrix of zeros
	Matrix<std::common_type_t<TA, TB>, M, P> result;
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

template<std::floating_point T, uint32_t M, uint32_t N, typename F>
Matrix<T, M, N> operator*(const Matrix<T, M, N>& matrix, const F& factor)
requires (std::is_arithmetic_v<F>)
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

template<std::floating_point T, uint32_t M, uint32_t N, typename F>
Matrix<T, M, N> operator*(const F& factor, const Matrix<T, M, N>& matrix)
requires (std::is_arithmetic_v<F>)
{
	return matrix * factor;
}

template<std::floating_point TA, std::floating_point TB, uint32_t M, uint32_t N>
Matrix<std::common_type_t<TA, TB>, M, N> operator+(const Matrix<TA, M, N>& mrxA, const Matrix<TB, M, N>& mrxB)
{
	Matrix<std::common_type_t<TA, TB>, M, N> result();
	for (unsigned rowIndex = 0; rowIndex < M; ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < N; ++columnIndex)
		{
			result[rowIndex][columnIndex] = mrxA[rowIndex][columnIndex] + mrxB[rowIndex][columnIndex];
		}
	}

	return result;
}

template<std::floating_point TA, std::floating_point TB, uint32_t M, uint32_t N>
Matrix<std::common_type_t<TA, TB>, M, N> operator-(const Matrix<TA, M, N>& mrxA, const Matrix<TB, M, N>& mrxB)
{
	Matrix<std::common_type_t<TA, TB>, M, N> result();
	for (unsigned rowIndex = 0; rowIndex < M; ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < N; ++columnIndex)
		{
			result[rowIndex][columnIndex] = mrxA[rowIndex][columnIndex] - mrxB[rowIndex][columnIndex];
		}
	}

	return result;
}

template<std::floating_point T, uint32_t M, uint32_t N>
Matrix<T, M, N> operator-(const Matrix<T, M, N>& mrxA)
{
	Matrix<T, M, N> result();
	for (unsigned rowIndex = 0; rowIndex < M; ++rowIndex)
	{
		for (unsigned columnIndex = 0; columnIndex < N; ++columnIndex)
		{
			result[rowIndex][columnIndex] = -mrxA[rowIndex][columnIndex];
		}
	}

	return result;
}

template<std::floating_point TA, std::floating_point TB, uint32_t M, uint32_t N>
bool operator==(const Matrix<TA, M, N>& mrxA, const Matrix<TB, M, N>& mrxB)
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

template<std::floating_point TA, std::floating_point TB, uint32_t M, uint32_t N>
bool operator!=(const Matrix<TA, M, N>& mrxA, const Matrix<TB, M, N>& mrxB)
{
	return !(mrxA == mrxB);
}
