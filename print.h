#include "matrix.h"

void PrintVectorFloat(Vector<float> v);
int GiveAbsInt(Matrix<int> A);
uint32_t Razryad(int max);
void PrintMatrixInt(Matrix<int> A);
float GiveMaxEl(Matrix<int> A);
float GiveMinEl(Matrix<int> A);
float Trace(Matrix<int> A);

int GiveAbsInt(Matrix<int> A) {
    int max = (A(0,0));
    int min = (A(0,0));
    for (uint32_t i=0; i<A.VertDim(); i++){
        for (uint32_t j=0; j<A.HorizDim(); j++){
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

uint32_t Razryad(int max) {
    uint32_t k =1; ///the number of razryadov
    if (max>=10) {while ((max/=10) > 0) k++;}
    return k;
}

void PrintMatrixInt(Matrix<int> A) {
    int max = GiveAbsInt(A);
    uint32_t k = Razryad(max),l = 1; ///the number of razryadov
    //while ((max/=10) > 0) k++;
    //std::cout << max<< " ";
    //std::cout << k<< std::endl;

    for (uint32_t i=0; i<A.VertDim(); i++){
        for (uint32_t j=0; j<A.HorizDim(); j++){
            int z = A(i,j);
            std::cout << z << " ";

            l=Razryad(z);
            while (k!=l){
                l++;
                std::cout << " ";
            }
        }
        std::cout<<std::endl;
    }

}

float GiveMaxEl(Matrix<int> A) {
    float max = (A(0,0));
    for (uint32_t i=0; i<A.VertDim(); i++){
        for (uint32_t j=0; j<A.HorizDim(); j++){
            if (A(i,j) > max) {
                max = A(i,j);
            }
        }}

    return max;
}

float GiveMinEl(Matrix<int> A) {
    int min = (A(0,0));
    for (uint32_t i=0; i<A.VertDim(); i++){
        for (uint32_t j=0; j<A.HorizDim(); j++){
            if (A(i,j) <min) {
                min = A(i,j);
            }
        }}
    return min;
}

float Trace(Matrix<int> A) {
    if (A.VertDim()!=A.HorizDim()){
        return -1;///here need to be exception, but it's not working...
    }
    else{
        float tr=0;
        for (uint32_t i=0; i<A.VertDim(); i++){
            tr+=A(i,i);
        }
        return tr;
    }}

void PrintVectorFloat(Vector<float> v) {
    for (uint32_t i=0; i<v.VertDim(); i++){
        std::cout << v(i)<< std::endl;
    }
}


