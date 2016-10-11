#ifndef FEMAREA_H
#define FEMAREA_H

#include "typedef.h"
#include <vector>
#include <list>

/*! struct FemArea - структура, описывающая МКЭ-область.
 *   границы по X: от x0 до x1
 *   границы по Y: от y0 до y1
 *   относительная магнитная проницаемость: muo
 *   значение плотности тока: j
 *   номер материала: nmat
 */
struct FemArea
{
    real x0, x1;    //!< границы по X: от x0 до x1.
    real y0, y1;    //!< границы по Y: от y0 до y1.
    int nmat;       //!< номер материала.
};

typedef std::vector<FemArea> FemAreaList;
typedef FemAreaList::iterator iFemAreaList;

#endif // FEMAREA_H
