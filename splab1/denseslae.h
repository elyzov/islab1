#ifndef DENSESLAE_H
#define DENSESLAE_H

#include "matrix.h"

class DenseSLAE
{
    Matrix<real> a;
    rvector f;
    rvector x;
    uint _dim;

    static constexpr real CMP_EPS = 1.e-14;
public:
    DenseSLAE(uint dim);
    real & A(uint i, uint j) {return a(i, j);}
    real & F(uint i) {return f.at(i);}
    real X(uint i) {return x.at(i);}
    void solve();
    void test();
};

class DenseSLAE_Error
{
    const char * str;
public:
    DenseSLAE_Error(const char * msg) : str(msg) {}
    void print_debug()
    {
        std::cerr << str << std::endl;
    }
};

#endif // DENSESLAE_H
