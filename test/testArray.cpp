#include <iostream>

#include "../src/array.h"


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
    std::cout << "Starting array test..." << std::endl;

    Array2f x(2,2,1.0);
    x(0,0) = 1.0;
    x(1,0) = 2.0;
    x(0,1) = 3.0;
    x(1,1) = 4.0;

    int numFailed = 0;
    numFailed += test(x(0,0) == 1.0, "operator(,) ");
    numFailed += test(x.pos(0,0).x == 0.5, "pos");

    numFailed += test(x.bilerp(0.0,0.0) == 1.0, "bilerp");
    numFailed += test(x.bilerp(2.0,0.0) == 2.0, "bilerp");    
    numFailed += test(x.bilerp(0.0,2.0) == 3.0, "bilerp");
    numFailed += test(x.bilerp(2.0,2.0) == 4.0, "bilerp");
    numFailed += test(x.bilerp(1.0,1.0) == (1+2+3+4)/4.0, "bilerp");
    numFailed += test(x.bilerp(1.0,0.5) == 1.5, "bilerp");
    
    numFailed += test(x.infNorm() == 4.0, "infNorm");
    
    x.reset();
    numFailed += test(x(0,0) == 0, "reset");
    
    Array2f y(2,2,1.0);
    y(0,0) = 1.0;
    y(1,0) = 2.0;
    y(0,1) = 3.0;
    y(1,1) = 4.0;

    x.swap(y);
    numFailed += test(x(0,0) == 1.0, "swap");

    y.copy(x);
    numFailed += test(x(0,0) == 1.0, "copy");

    numFailed += test(x.dot(y) == 30, "dot");

    x.add(y);
    numFailed += test(x(0,0) == 2.0, "add");

    x.scaleAndAdd(2.0,y);
    numFailed += test(x(0,0) == 5.0, "scaleAndAdd");

    x.resize(3,3);
    numFailed += test(x.nx() == 3 && x.ny() == 3, "resize");    
    
    FaceArray2Xf u(2.0,2.0,1.0);
    FaceArray2Yf v(2.0,2.0,1.0);

    u.face<LEFT>(0,0) = 3.0;
    u.face<LEFT>(1,0) = 4.0;
    u.face<RIGHT>(1,0) = 5.0;
    u.face<LEFT>(0,1) = 8.0;
    u.face<LEFT>(1,1) = 9.0;
    u.face<RIGHT>(1,1) = 10.0;

    v.face<BOTTOM>(0,0) = 1.0;
    v.face<BOTTOM>(1,0) = 2.0;
    v.face<BOTTOM>(0,1) = 6.0;
    v.face<BOTTOM>(1,1) = 7.0;
    v.face<TOP>(0,1) = 11.0;
    v.face<TOP>(1,1) = 12.0;

    numFailed += test(u.pos<LEFT>(0,0).x == 0.0, "pos face");
    numFailed += test(u.pos<RIGHT>(1,1).x == 2.0, "pos face");
    numFailed += test(u.pos<RIGHT>(1,1).y == 1.5, "pos face");
    
    numFailed += test(u.center(0,0) == 7.0/2.0, "bilerp center");
    numFailed += test(u.bilerp(0.0,0.0) == 3.0, "bilerp face");
    numFailed += test(u.bilerp(1.0,0.0) == 4.0, "bilerp face");
    numFailed += test(u.bilerp(2.0,0.5) == 5.0, "bilerp face");
    numFailed += test(u.bilerp(0.0,2.0) == 8.0, "bilerp face");
    numFailed += test(u.bilerp(2.0,2.0) == 10.0, "bilerp face");
    numFailed += test(u.bilerp(1.0,1.0) == 13/2.0, "bilerp face");
    numFailed += test(u.bilerp(0.5,1.0) == (3+4+8+9)/4.0, "bilerp face");
    numFailed += test(u.bilerp(1.0,1.5) == (6+7+11+12)/4.0, "bilerp face");
    numFailed += test(u.bilerp(1.5,2.0) == 19/2.0, "bilerp face");

    numFailed += test(v.pos<BOTTOM>(0,0).x == 0.5, "pos face");
    numFailed += test(v.pos<BOTTOM>(1,1).x == 1.5, "pos face");
    numFailed += test(v.pos<TOP>(1,1).y == 2.0, "pos face");
    
    numFailed += test(v.center(0,0) == 7.0/2.0, "bilerp center");
    numFailed += test(v.bilerp(0.0,0.0) == 1.0, "bilerp face");
    numFailed += test(v.bilerp(1.0,0.0) == 1.5, "bilerp face");
    numFailed += test(v.bilerp(2.0,0.5) == 9.0/2.0, "bilerp face");
    numFailed += test(v.bilerp(0.0,2.0) == 11.0, "bilerp face");
    numFailed += test(v.bilerp(2.0,2.0) == 12.0, "bilerp face");
    numFailed += test(v.bilerp(1.0,1.0) == 13/2.0, "bilerp face");
    numFailed += test(v.bilerp(0.5,1.0) == 6.0, "bilerp face");
    numFailed += test(v.bilerp(1.0,1.5) == 9.0, "bilerp face");
    numFailed += test(v.bilerp(1.5,2.0) == 12.0, "bilerp face");
    numFailed += test(v.bilerp(1.0,0.5) == 4.0, "bilerp face");

    CornerArray2f c(2,2,1.0);
    c.corner<BOTTOM_LEFT>(0,0) = 1.0;
    c.corner<BOTTOM_LEFT>(1,0) = 2.0;
    c.corner<BOTTOM_RIGHT>(1,0) = 3.0;
    c.corner<BOTTOM_LEFT>(0,1) = 4.0;
    c.corner<BOTTOM_LEFT>(1,1) = 5.0;
    c.corner<BOTTOM_RIGHT>(1,1) = 6.0;
    c.corner<TOP_LEFT>(0,1) = 7.0;
    c.corner<TOP_LEFT>(1,1) = 8.0;
    c.corner<TOP_RIGHT>(1,1) = 9.0;

    
    numFailed += test(c.pos<BOTTOM_LEFT>(0,0).x == 0.0, "pos face");
    numFailed += test(c.pos<BOTTOM_RIGHT>(1,1).x == 2.0, "pos face");
    numFailed += test(c.pos<TOP_RIGHT>(1,1).y == 2.0, "pos face");

    numFailed += test(c.bilerp(0.0,0.0) == 1.0, "bilerp face");
    numFailed += test(c.bilerp(2.0,0.0) == 3.0, "bilerp face");
    numFailed += test(c.bilerp(1.0,0.0) == 2.0, "bilerp face");
    numFailed += test(c.bilerp(1.0,1.0) == 5.0, "bilerp face");

    numFailed += test(c.bilerp(0.0,2.0) == 7.0, "bilerp face");
    numFailed += test(c.bilerp(2.0,2.0) == 9.0, "bilerp face");
    numFailed += test(c.bilerp(2.0,1.5) == (9+6)/2.0, "bilerp face");
    numFailed += test(c.bilerp(1.5,1.5) == (5+6+8+9)/4.0, "bilerp face");
    numFailed += test(c.center(1,1) == (5+6+8+9)/4.0, "center");
    
    std::cout << "Number of failed tests: " << numFailed << std::endl;    
}
