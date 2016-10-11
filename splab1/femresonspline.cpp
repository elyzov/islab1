#include "femresonspline.h"
#include "femerror.h"

FemResOnSpline::FemResOnSpline(const std::string & value, const std::string & points)
{
    if (spline.build(value) < 0) throw Fem_OpenFile_Error(value);;
    inputPoints(points);
}

void FemResOnSpline::inputPoints(const std::string & fname)
{
    std::ifstream fin(fname.c_str());
    if (!fin) throw Fem_OpenFile_Error(fname);
    fin.exceptions(std::ifstream::failbit |
                   std::ifstream::badbit  |
                   std::ifstream::eofbit);

    uint num;
    fin >> num;
    xySources.resize(num);
    powSources.resize(num);
    for (uint i = 0; i < num; ++i)
    {
        fin >> xySources.at(i).first.x >> xySources.at(i).first.y;
        fin >> xySources.at(i).second.x >> xySources.at(i).second.y;
        fin >> powSources.at(i);
    }

    fin >> num;
    xyReceivers.resize(num);
    srcReceivers.resize(num);
    for (uint i = 0; i < num; ++i)
    {
        fin >> xyReceivers.at(i).first.x
            >> xyReceivers.at(i).first.y;
        fin >> xyReceivers.at(i).second.x
            >> xyReceivers.at(i).second.y;
        fin >> srcReceivers.at(i);
    }
}

int FemResOnSpline::diffSolInPointsPair(rvector & sol)
{
    if (sol.empty()) sol.resize(xyReceivers.size());
    for (uint i = 0, iend = xyReceivers.size(); i != iend; ++i)
    {
        FemPoint & pM = xyReceivers.at(i).first;
        FemPoint & pN = xyReceivers.at(i).second;
        FemPoint & pA = xySources.at(srcReceivers.at(i)).first;
        FemPoint & pB = xySources.at(srcReceivers.at(i)).second;
        sol.at(i) = spline.value(len(pA, pN)) - spline.value(len(pA, pM)) -
                    spline.value(len(pB, pN)) + spline.value(len(pB, pM));
        sol.at(i) *= powSources.at(srcReceivers.at(i));
    }
    return 0;
}
