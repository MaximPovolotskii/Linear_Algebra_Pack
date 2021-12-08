#ifndef LINALG_VECT_H
#define LINALG_VECT_H

#include <vector>
#include <exception>
#include "matrix.h"

int* PowersTwo(int a){///4 is equal to [0,0,1]
    int a_2[a];
    int i = 0;
    while (a > 0) {
        a_2[i] = a % 2;
        a = a / 2;
        i++;
    }
    std::cout<<a_2;
}
#endif //LINALG_VECT_H
