#ifndef FEMRESONSPLINE_H
#define FEMRESONSPLINE_H

#include "typedef.h"
#include "cubicspline.h"
#include "fempoint.h"

class FemResOnSpline
{
    CubicSpline spline;
    FemPointPairList xySources;
    rvector powSources;
    FemPointPairList xyReceivers;
    intVector srcReceivers;


public:
    FemResOnSpline(const std::string & value, const std::string & points);
    int diffSolInPointsPair(rvector & sol);

private:
    void inputPoints(const std::string &fname);
};

#endif // FEMRESONSPLINE_H
