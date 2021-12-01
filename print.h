#include "matrix.h"
int give_abs_int(Matrix<int> A, int n, int m){
    int max = (A(0,0));
    int min = (A(0,0));
    for (uint32_t i=0; i<n; i++){
        for (uint32_t j=0; j<m; j++){
            if (A(i,j) <min) {
                min = A(i,j);
            }
                if (A(i,j) > max) {
                    max = A(i,j);
                }
    }}
    if (abs(min)>abs(max)){
    return abs(min);}
    else {

        return abs(max);}
}

uint32_t razryad(int max){
    uint32_t k =1; ///the number of razryadov
    if (max>=10) {while ((max/=10) > 0) k++;}
    return k;
}

void print_matrix_int(Matrix<int> A, int n, int m){
    int max = give_abs_int(A, n, m);
    uint32_t k = razryad(max),l = 1; ///the number of razryadov
    //while ((max/=10) > 0) k++;
    //std::cout << max<< " ";
    //std::cout << k<< std::endl;

    for (uint32_t i=0; i<n; i++){
        for (uint32_t j=0; j<m; j++){
            int z = A(i,j);
            std::cout << z << " ";

            l=razryad(z);
            while (k!=l){
                l++;
                std::cout << " ";
            }
        }
        std::cout<<std::endl;
    }

}
float give_max_el(Matrix<int> A, int n, int m){
    float max = (A(0,0));
    for (uint32_t i=0; i<n; i++){
        for (uint32_t j=0; j<m; j++){
            if (A(i,j) > max) {
                max = A(i,j);
            }
        }}

        return max;
}

float give_min_el(Matrix<int> A, int n, int m){
    int min = (A(0,0));
    for (uint32_t i=0; i<n; i++){
        for (uint32_t j=0; j<m; j++){
            if (A(i,j) <min) {
                min = A(i,j);
            }
        }}
        return min;
}

float give_trace(Matrix<int> A, int n, int m){
    if (n!=m){
        return -1;///here need to be exception, but it's not working...
    }
    else{
        float tr=0;
        for (uint32_t i=0; i<n; i++){
            tr+=A(i,i);
        }
        return tr;
}}