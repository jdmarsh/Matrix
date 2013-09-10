#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix {
public:
    Matrix(unsigned, unsigned);
    unsigned cols();
    unsigned rows();
    float determinant();
    float trace();
    Matrix inverse();
    Matrix transpose();
    Matrix REF();
    Matrix RREF();
    unsigned rank();
    void scale(float);

    std::vector<float>& operator[](unsigned);
    void operator+=(Matrix);
    void operator-=(Matrix);
private:
    std::vector<std::vector<float>> matrix;

    void rowdivide(unsigned, float);
    void rowswap(unsigned, unsigned);
    void rowminus(unsigned, unsigned, float);
    Matrix newmrx(unsigned, unsigned);
};

Matrix id(unsigned);
Matrix operator*(Matrix, Matrix);
Matrix operator+(Matrix, Matrix);
Matrix operator-(Matrix, Matrix);
Matrix operator-(Matrix);
bool operator==(Matrix, Matrix);
bool operator!=(Matrix, Matrix);

#endif