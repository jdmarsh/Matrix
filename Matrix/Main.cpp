#include "Main.h"

int main() {
    Matrix A(3, 2);
    A[0][0] = 6;
    A[0][1] = 2;
    A[1][0] = 0;
    A[1][1] = 4;
    A[2][0] = 5;
    A[2][1] = 10;
    Matrix E = A.RREF();
    Matrix D = -A;
    std::cout << D[0][0] << ",";
    std::cout << D[0][1] << std::endl;
    std::cout << D[1][0] << ",";
    std::cout << D[1][1] << std::endl;
    std::cout << D[2][0] << ",";
    std::cout << D[2][1];
    return 0;
}