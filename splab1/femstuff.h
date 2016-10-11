#ifndef FEMSTUFF_H
#define FEMSTUFF_H

#include "typedef.h"
#include "cubicspline.h"
#include <map>

/*! struct FemMat - описание материала для МКЭ-задачи.
 *   относительная магнитная проницаемость: muo[i]
 *   значение плотности тока: j[i]
 *   номер материала: nmat[i]
 */
class FemStuff
{
//    bool linear_muo;
    real _sig;                  //!< проводимость.
//    CubicSpline _muo_spline;    //!< для нелинейного случая.

public: /*! Интерфейс. */
    FemStuff(real sig = 0) : _sig(sig)
    {
//        linear_muo = true;
    }
//    void test_nonlinearty(const std::string & file)
//            {if (_muo_spline.build(file) == 0) linear_muo = false;}

/*! Get-методы. */
//    real sig(real B, ubyte der = 0) {return _muo_spline.value(B, der);}
    real sig() {return _sig;}
//    bool is_linear() {return linear_muo;}
//    real greatest_muo_arg() {return _muo_spline.right_boundary();}

/*! Set-методы. */
    void set_sig(real sig) {_sig = sig;}
};

typedef std::map<int, FemStuff> FemStuffList;
typedef FemStuffList::iterator iFemStuffList;
typedef FemStuffList::const_iterator ciFemStuffList;

#endif // FEMSTUFF_H
