#include <iostream>

#include "../src/sparse.h"

bool printPassed = false;

int test(bool cond, const char * msg)
{
    if (cond) {
        if (printPassed) {
            std::cout << msg << " ... PASSED" << std::endl;
        }
        return 0;
    } else {
        std::cout << msg << " ... FAILED" << std::endl;
        return 1;
    }
}

int main(int argc, char *argv[]) {
    std::cout << "Starting sparse matrix test..." << std::endl;

    SparseLaplacianMatrix<float> A(3,3);

    A.value<CENTER>(1,1) = 1.0;
    A.value<LEFT>(1,1) = 1.0;
    A.value<RIGHT>(1,1) = 2.0;
    A.value<BOTTOM>(1,1) = 3.0;
    A.value<TOP>(1,1) = 4.0;


    A.value<CENTER>(0,0) = 1.0;
    A.value<TOP>(0,0) = 2.0;
    A.value<RIGHT>(0,0) = 3.0;

    Array2f x(3,3,1.0);
    x.set(1.0);

    int numFailed = 0;

    numFailed += test(A.mult(x,0,0) == 6.0, "mult");
    numFailed += test(A.mult(x,1,1) == 11.0, "mult");
}
