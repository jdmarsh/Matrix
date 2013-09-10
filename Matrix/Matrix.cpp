#include "Matrix.h"

Matrix::Matrix(unsigned i, unsigned j) {
    matrix.resize(i);
    for (auto& it : matrix) {
        it.resize(j);
    }
}

unsigned Matrix::cols() {
    if (matrix.size() == 0) {
        return 0;
    }
    return matrix[0].size();
}

unsigned Matrix::rows() {
    return matrix.size();
}

float Matrix::determinant() {
    float determinant = 0;
    if (cols() != rows() || rows() == 0) {
        return 0;
    }
    if (rows() == 1) {
        return (matrix[0][0]);
    }

    for (unsigned j = 0; j < rows(); ++j) {
        float a = matrix[0][j] * newmrx(0, j).determinant();
        determinant += ((j % 2) ? -a : a);
    }
    return determinant;
}

float Matrix::trace() {
    float tr = 0;
    if (rows() == cols()) {
        for (unsigned i = 0; i < rows(); ++i) {
            tr += matrix[i][i];
        }
    }
    return tr;
}

Matrix Matrix::inverse() {
    float det = determinant();
    if (det == 0) {
        return id(0);
    }
    if (matrix.size() == 1) {
        Matrix X = id(1);
        X.scale(1.0f / matrix[0][0]);
        return X;
    }
    Matrix A(rows(), cols());
    for (unsigned i = 0; i < rows(); ++i) {
        for (unsigned j = 0; j < cols(); ++j) {
            float a = newmrx(i, j).determinant();
            A.matrix[i][j] = (((i + j) % 2) ? -a : a);
        }
    }
    A.transpose().scale(1.0f / det);
    return A;
}

Matrix Matrix::transpose() {
    Matrix B(cols(), rows());
    for (unsigned i = 0; i < cols(); ++i) {
        for (unsigned j = 0; j < rows(); ++j) {
            B.matrix[i][j] = matrix[j][i];
        }
    }
    return B;
}

Matrix Matrix::REF() {
    Matrix M(rows(), cols());
    M.matrix = matrix;
    unsigned x = 0;
    for (unsigned j = 0; j < cols(); ++j) {
        for (unsigned i = x; i < rows(); ++i) {
            if (M.matrix[i][j]) {
                M.rowswap(i, x);
                M.rowdivide(x, M.matrix[x][j]);
                for (unsigned a = x + 1; a < rows(); ++a) {
                    M.rowminus(a, x, M.matrix[a][j]);
                }
                x++;
                break;
            }
        }
    }
    return M;
}

Matrix Matrix::RREF() {
    Matrix M = REF();
    for (unsigned i = 0; i < M.rows(); ++i) {
        for (unsigned j = 0; j < M.cols(); ++j) {
            if (M.matrix[i][j]) {
                for (unsigned a = 0; a < i; ++a) {
                    M.rowminus(a, i, M.matrix[a][j]);
                }
                break;
            }
        }
    }
    return M;
}

unsigned Matrix::rank() {
    Matrix M = RREF();
    unsigned x = 0;
    for (unsigned i = 0; i < M.rows(); ++i) {
        for (unsigned j = 0; j < M.cols(); ++j) {
            if (M.matrix[i][j]) {
                x++;
                break;
            }
        }
    }
    return x;
}

void Matrix::scale(float x) {
    for (unsigned i = 0; i < rows(); ++i) {
        for (unsigned j = 0; j < cols(); ++j) {
            matrix[i][j] *= x;
        }
    }
}

std::vector<float>& Matrix::operator[](unsigned i){
    return matrix[i];
}

void Matrix::operator+=(Matrix A) {
    if ((rows() == A.rows()) && (cols() == A.cols())) {
        for (unsigned i = 0; i < rows(); ++i) {
            for (unsigned j = 0; j < cols(); ++j) {
                matrix[i][j] += A[i][j];
            }
        }
    }
}

void Matrix::operator-=(Matrix A) {
    if ((rows() == A.rows()) && (cols() == A.cols())) {
        for (unsigned i = 0; i < rows(); ++i) {
            for (unsigned j = 0; j < cols(); ++j) {
                matrix[i][j] -= A[i][j];
            }
        }
    }
}

void Matrix::rowdivide(unsigned i, float x) {
    for (unsigned a = 0; a < cols(); ++a) {
        matrix[i][a] /= x;
    }
}

void Matrix::rowswap(unsigned i, unsigned j) {
    std::vector<float> vec;
    vec = matrix[i];
    matrix[i] = matrix[j];
    matrix[j] = vec;
}

void Matrix::rowminus(unsigned i, unsigned j, float x) {
    for (unsigned a = 0; a < cols(); ++a) {
        matrix[i][a] -= (matrix[j][a] * x);
    }
}

Matrix Matrix::newmrx(unsigned i, unsigned j) {
    if (rows() == 0) return id(0);
    Matrix mrx(rows(), cols());
    mrx.matrix = matrix;
    for (unsigned a = 0; a < rows(); ++a) {
        mrx.matrix[a].erase(mrx.matrix[a].begin() + j);
    }
    mrx.matrix.erase(mrx.matrix.begin() + i);
    return mrx;
}

Matrix id(unsigned dim) {
    Matrix i(dim, dim);
    for (unsigned a = 0; a < dim; ++a) {
        i[a][a] = 1;
    }
    return i;
}

Matrix operator*(Matrix A, Matrix B) {
    if (A.cols() != B.rows()) {
        return id(0);
    }
    Matrix C(A.rows(), B.cols());
    for (unsigned i = 0; i < C.rows(); ++i) {
        for (unsigned j = 0; j < C.cols(); ++j) {
            float sum = 0;
            for (unsigned k = 0; k < A.cols(); ++k) {
                sum += (A[i][k] * B[k][j]);
            }
            C[i][j] = sum;
        }
    }
    return C;
}

Matrix operator+(Matrix A, Matrix B) {
    if ((A.rows() != B.rows()) || (A.cols() != B.cols())) {
        return id(0);
    }
    for (unsigned i = 0; i < A.rows(); ++i) {
        for (unsigned j = 0; j < A.cols(); ++j) {
            A[i][j] += B[i][j];
        }
    }
    return A;
}

Matrix operator-(Matrix A, Matrix B) {
    if ((A.rows() != B.rows()) || (A.cols() != B.cols())) {
        return id(0);
    }
    for (unsigned i = 0; i < A.rows(); ++i) {
        for (unsigned j = 0; j < A.cols(); ++j) {
            A[i][j] -= B[i][j];
        }
    }
    return A;
}

Matrix operator-(Matrix A) {
    for (unsigned i = 0; i < A.rows(); ++i) {
        for (unsigned j = 0; j < A.cols(); ++j) {
            A[i][j] = -A[i][j];
        }
    }
    return A;
}

bool operator==(Matrix A, Matrix B) {
    if ((A.rows() != B.rows()) || (A.cols() != B.cols())) {
        return 0;
    }
    for (unsigned i = 0; i < A.rows(); ++i) {
        for (unsigned j = 0; j < A.cols(); ++j) {
            if (A[i][j] != B[i][j]) {
                return 0;
            }
        }
    }
    return 1;
}

bool operator!=(Matrix A, Matrix B) {
    return !(A == B);
}