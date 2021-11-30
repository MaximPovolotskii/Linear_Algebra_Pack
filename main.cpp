#include <iostream>
#include "matrix.h"

int main() {
    int c[4] = {1, 2, 3, 4};
    Matrix<int> a(c, 2, 2);
    Matrix<int> b = a.SmartMult(a);
    std::cout << b(0, 0) << " " << b(0, 1) << " " << b(1, 0) << " " << b(1, 1) << std::endl;
    return 0;
}
