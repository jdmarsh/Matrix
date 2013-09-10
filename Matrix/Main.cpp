#include "Main.h"

int main() {
    Matrix M{ { 1, 2 }, { 3, 4 }, { 5, 6 } };
    Matrix A(3, 2);
    A[0][0] = 6;
    A[0][1] = 2;
    A[1][0] = 0;
    A[1][1] = 4;
    A[2][0] = 5;
    A[2][1] = 10;
    Matrix E = A.RREF();
    Matrix D = -A;
    std::cout << M[0][0] << ",";
    std::cout << M[0][1] << std::endl;
    std::cout << M[1][0] << ",";
    std::cout << M[1][1] << std::endl;
    std::cout << M[2][0] << ",";
    std::cout << M[2][1];
    return 0;
}