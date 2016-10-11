#ifndef FEMGRID_H
#define FEMGRID_H

#include "typedef.h"
#include <vector>

/*! string FemGrid - структура описывающая МКЭ-сетку.
 *  сoord[i] - граница,
 *  hMin[i] - минимальный шаг на участке от сoord[i-1] до coord[i],
 *  dh[i] - коэффициент разрядки на участке от coord[i-1] до coord[i],
 *  sh[i] - знак коэффициент разрядки на участке от coord[i-1] до coord[i]:
             если -1, то шаги увеличиваются сверху вниз,
             если  1, то шаги увеличиваются снизу вверх.
 */
struct FemGrid
{
    real coord; //!< граница.
    real hMin;  //!< минимальный шаг.
    real dh;    //!< коэффициент разрядки.
    real sh;    //!< коэффициент разрядки.

    FemGrid & operator *= (real a)
    {
        coord *= a;
        hMin *= a;
        return *this;
    }
};

typedef std::vector<FemGrid> FemGridList;
typedef FemGridList::iterator iFemGridList;
typedef FemGridList::const_iterator ciFemGridList;

#endif // FEMGRID_H
