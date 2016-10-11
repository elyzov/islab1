#ifndef FEMPOINT_H
#define FEMPOINT_H

#include "typedef.h"
#include <vector>
#include <iostream>

struct FemPoint
{
    real x, y;      //!< координаты точки.
    FemPoint(real coordX = 0, real coordY = 0);
    void print(){std::cout << "(" << x << "," << y << ")" << std::endl;}
};

real len(const FemPoint & a, const FemPoint & b);
real len(const FemPoint & a);
FemPoint mid(const FemPoint & a, const FemPoint & b);

typedef std::vector<FemPoint> FemPointList;
typedef FemPointList::iterator iFemPointList;
typedef std::vector<std::pair<FemPoint, real> > pointRealVector;
typedef std::pair<FemPoint, FemPoint> FemPointPair;
typedef std::vector<FemPointPair> FemPointPairList;
typedef FemPointPairList::iterator iFemPointPairList;
typedef std::vector<FemPointPairList> FemPointPairList2d;

#endif // FEMPOINT_H
