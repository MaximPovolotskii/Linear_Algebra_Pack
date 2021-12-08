#include <iostream>
#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch_amalgamated.hpp"
#include "matrix.h"
#include "vector.h"
#include <cstdlib>
#include <cmath>
#include "linear_system_solve.h"

/*
TEST_CASE("operator+, operator=, .Shape() without random numbers") {
    int n=3, m=4, N = n*m;
    int max;
    int a[12] {1000, 2530, 0, 0, 2, 260, 3, 4, 70, 0, -8, 1};
    int b[12] {1, -600, 265, 305, 2, 0, 3, 4, 70, 0, -8, 1};
    int c[12] {1001, 1930, 265, 305, 4, 260, 6, 8, 140, 0, -16, 2};

    Matrix<int> A(a, n, m);
    Matrix<int> B(b, n, m);
    Matrix<int> tested = A+B;
    Matrix<int> C(c, n, m);
    Matrix<int> C2 = B+A;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            REQUIRE(tested(i,j) == C(i,j));
            REQUIRE(C2(i,j) == C(i,j));
        }
    }
    REQUIRE(tested.Shape() == A.Shape());
    REQUIRE(tested.Shape() == B.Shape());
}

TEST_CASE("operator+, operator=, .Shape(), operator-, operator() with random elements in range(0,32767)") {
    int n=round(rand() * (20) /3000)+1, m=round(rand() * (20) /3000)+1;
    std::cout<<n<<m;
    int a[m*n];
    int b[n*m];
    int c[n*m];
    int d[n*m];
    for (int i=0;i<n*m;i++){
        a[i] = rand()*(-1)^rand();
        b[i] = rand()*(-1)^rand();
        c[i] = a[i]+b[i];
        d[i] = a[i]-b[i];
    }

    Matrix<int> A(a, n, m);
    Matrix<int> B(b, n, m);
    Matrix<int> tested = A+B;
    Matrix<int> C(c, n, m);
    Matrix<int> C2 = B+A;
    Matrix<int> testD = A-B;
    Matrix<int> D(d, n, m);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            REQUIRE(tested(i,j) == C(i,j));
            REQUIRE(C2(i,j) == C(i,j));
            REQUIRE(testD(i,j) == D(i,j));
        }
    }
    REQUIRE(tested.Shape() == A.Shape());
    REQUIRE(tested.Shape() == B.Shape());
}

TEST_CASE("operator* without random elements") {
    int n=3, m=5, k = 2;
    int a[15]{7,92,3,0,4,7,5,4,6,9,1,2,7,6,1};
    int b[5*2]{43,56,0,5,0,61,0,2,0,1};
    int c[3*2]{301,1039,301,682,43,506};

    Matrix<int> A(a, n, m);
    Matrix<int> B(b, m, k);
    Matrix<int> tested = A*B;
    Matrix<int> C(c, n, k);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            REQUIRE(tested(i,j) == C(i,j));

        }
    }
    REQUIRE(tested.VertDim() == A.VertDim());
    REQUIRE(tested.HorizDim() == B.HorizDim());

}


TEST_CASE("operator==,operator!= with random elements in range(0,32767)") {
    int n=round(rand() * (20) /3000)+1, m=round(rand() * (20) /3000)+1;
    int a[m*n];
    int b[n*m];

    for (int i=0;i<n*m;i++){
        a[i] = rand()*(-1)^rand();
        b[i] = rand()*(-1)^rand();
        if (a[i]==b[i]) b[i]-=1;

    }

    Matrix<int> A(a, n, m);
    Matrix<int> B(b, n, m);
    Matrix<int> tested = A+B;
    Matrix<int> C(a, n, m);


    REQUIRE(A == C);
    REQUIRE(A != B);
    REQUIRE(A==A);
    REQUIRE(B!=A);
    REQUIRE(B==B);
}


TEST_CASE("matrix in power with random elements") {
    int n=4;
    int a2[16];
    int a4[16];
    for (int i=0;i<n*n;i++){
        a2[i] = rand()*(-1)^rand();
        a4[i] = rand()*(-1)^rand();
    }
    int a[16]{4,56,0,32,2,56,21,0,58,9,6,52,5,1,2,5};

    Matrix<int> A(a, n, n);
    Matrix<int> D(a2, n, n);
    int c[9] {4,3,2, 0,7,5, 0, 0,1};
    Matrix<int> Id({1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1}, n, n);
    Matrix<int> C(c, 3, 3);

    REQUIRE(A.Pow(A,2) == A*A);
    REQUIRE(A*A*A*A == A.Pow(A,4));
    REQUIRE(C*C*C == C.Pow(C,3));
    REQUIRE(C*C*C*C*C == C.Pow(C,5));
    REQUIRE(C*C*C*C*C*C == C.Pow(C,6));

    REQUIRE(D.Pow(D,2) == D*D);
    REQUIRE(D*D*D*D == D.Pow(D,4));
    REQUIRE(D*D*D == D.Pow(D,3));
    REQUIRE(D*D*D*D*D == D.Pow(D,5));
    REQUIRE(D*D*D*D*D*D == D.Pow(D,6));
    REQUIRE(D.Pow(D,1) == D);
}*/
TEST_CASE("determinant and rank") {
    float n=5;
    float a2[25]{54,	8,	9,	5,	2,
    2,	5,	55,	5,	5,
    6,	5,	5,	1,	2,
    3,	3,	5,	4,	4,
    52,	2,	2,	0,	0};


    float a[16]{4,56,0,32,2,56,21,0,58,9,6,52,5,1,2,5};

    Matrix<float> A(a, 4, 4);
    Matrix<float> D(a2, n, n);
    float c[9] {4,3,2, 0,7,5, 0, 0,1};
    Matrix<float> Id({1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1}, 4, 4);
    Matrix<float> C(c, 3, 3);
    std::cout<<Id.Determinant();

    float k[9]{2,  12,  0, 1,  21,  2, 1,  0,  2};
    Matrix<float> K(k, 3, 3);

    REQUIRE(std::abs(K.Determinant() - 84) < 0.01);
    REQUIRE(std::abs(D.Determinant() + 143702) < 0.3);
    REQUIRE(std::abs(C.Determinant() - 28 ) < 0.01);
    REQUIRE(std::abs(Id.Determinant() - 1 ) < 0.01);

    float m[16]{1,1,  1,  1,
              1,  1,  11,  1,
              1,  1,  1,  1,
              1,  1,  1,  1};
    Matrix<float> M(m, 4, 4);
    REQUIRE(M.Determinant() == 0);

    REQUIRE(K.Rank() == 3);
    REQUIRE(D.Rank() == 5);
    REQUIRE(C.Rank() == 3);
    REQUIRE(Id.Rank() == 4);
    REQUIRE(M.Rank() == 2);

    float nf[16]{1, 2,  1,  1,
                1,  2,  1,  1,
                1,  2,  1,  1,
                1,  2,  1,  1};
    Matrix<float> N(nf, 4, 4);
    REQUIRE(N.Rank() == 1);
    Matrix<float> NUll(4, 4, NULLMATRIX);
    REQUIRE(NUll.Rank() == 0);

}
TEST_CASE("normalize") {
    float v1[3]{1,1,1};
    float v2[3]{1,2,3};
    float v3[5]{1,0,1,0,1};
    float v4[9] {0,0,1,0,0,0,0,0,1};
    Vector<float> V1(v1,3);
    Vector<float> V2(v2,3);
    Vector<float> V3(v3,5);
    Vector<float> V4(v4,9);
    V1 = Normalize(V1);
    V2 = Normalize(V2);
    V3 = Normalize(V3);
    V4 = Normalize(V4);
    REQUIRE(std::abs(V1(1)*V1(1)+V1(2)*V1(2)+V1(0)*V1(0)-1) <0.2);
    REQUIRE(std::abs(V2(1)*V2(1)+V2(2)*V2(2)+V2(0)*V2(0)-1) <0.2);
    REQUIRE(std::abs(V3(1)*V3(1)+V3(2)*V3(2)+V3(0)*V3(0)+V3(3)*V3(3)+V3(4)*V3(4)-1) <0.2);
}

TEST_CASE("TRANSPOSE") {

    float a2[25]{54,  8,  9,  5,  2,
                 2,  5,  55,  5,  5,
                 6,  5,  5,  1,  2,
                 3,  3,  5,  4,  4,
                 52,  2,  2,  0,  0};


    float a[16]{4,56,0,32,2,56,21,0,58,9,6,52,5,1,2,5};

    Matrix<float> A(a, 4, 4);
    Matrix<float> Dt(a2, 5, 5);
    Matrix<float> D=Dt;
    Dt.Transpose();
    float c[9] {4,3,2, 0,7,5, 0, 0,1};
    Matrix<float> Id({1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1}, 4, 4);
    Matrix<float> Ct(c, 3, 3);
    Matrix<float> C=Ct;
    Ct.Transpose();
    std::cout<<Id.Determinant();

    float k[9]{2,  12,  0, 1,  21,  2, 1,  0,  2};
    Matrix<float> Kt(k, 3, 3);
    Matrix<float> K=Kt;
    Kt.Transpose();

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            REQUIRE(K(i,j) == Kt(j,i));
            REQUIRE(Ct(i,j) == C(j,i));
        }
    }

    int n=3, m=4, N = 12;
    float f[12] {1000, 2530, 0, 0, 2, 260, 3, 4, 70, 0, 8, 1};
    Matrix<float> Ft(f, n, m);
    Matrix<float> F=Ft;
    Ft.Transpose();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            REQUIRE(F(i,j) == Ft(j,i));
        }
    }
}

TEST_CASE("Linear System solve") {
    float a1[25]{54,  8,  9,  5,  2,
                 2,  5,  55,  5,  5,
                 6,  5,  5,  1,  2,
                 3,  3,  5,  4,  4,
                 52,  2,  2,  0,  0};
    float v1[5]{54, 5, 5, 4, 0};

    Matrix<float> D1(a1, 5, 5);
    Vector<float> V1(v1, 5);

    auto p = LinearSystemSolve(D1, V1);

}
