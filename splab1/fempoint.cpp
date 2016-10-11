#include "fempoint.h"
#include <math.h>

FemPoint::FemPoint(real coordX, real coordY)
    : x(coordX), y(coordY)
{
}

real len(const FemPoint & a, const FemPoint & b)
{
    return sqrt(pow(b.x - a.x, 2) + pow(b.y - a.y, 2));
}

FemPoint mid(const FemPoint & a, const FemPoint & b)
{
    FemPoint m;
    m.x = (a.x + b.x) / 2;
    m.y = (a.y + b.y) / 2;
    return m;
}

