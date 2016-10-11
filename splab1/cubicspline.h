#ifndef CUBICSPLINE_H
#define CUBICSPLINE_H

#include "typedef.h"
#include <fstream>
#include <string>
#include "cubicsplineerror.h"
#include <limits>
#include <stdio.h>

/*! class CubicSpline - кубический сплайн дефекта 1. */
class CubicSpline
{
    struct tuple
    {
        real a, b, c, d, x, f;
        bool operator < (const tuple & op) const
        {
            return x < op.x;
        }
    };

    rvector f;
    std::vector<tuple> spline;
    rvector2d di;
    int n;
public: /*! Интерфейс. */
    CubicSpline();
    void clear();
    real value(real x, ubyte derivative = 0);
    int build(const std::string & file);
    real left_boundary();
    real right_boundary();

private:
    void sol_di_slau();
    void fill_matr();
    void calc_kof();
};

#endif // CUBICSPLINE_H
