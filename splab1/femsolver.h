#ifndef FEMSOLVER_H
#define FEMSOLVER_H

#include <iostream>
#include <fstream>
#include <string>
#include "femnamespace.h"
#include "typedef.h"
#include "femerror.h"
#include "fempoint.h"
#include "femstuff.h"
#include "femdatgen.h"
#include "slau.h"
#include "cubicspline.h"

using namespace FEM;

/*! class FemSolver - решатель уравнения МКЭ. */
class FemSolver
{
private:
/*! -------------------------------------------------
 *                  Входные данные                  |
 *  -------------------------------------------------
 */
    int kuzlov;         //!< число узлов сетки.
    int ktr;            //!< число конечных элементов.
    int kt1;            //!< число узлов с первыми краевыми условиями.

    intVector2d nvtr;     //!< массив конечных элементов (с номером материала).
    FemPointList xy;      //!< массив координат.
    intVector nvtrFst;    //!< массив первых краевых.
    FemStuffList stuff;   //!< каталог материалов.

/*! -------------------------------------------------
 *              Вспомагательные данные              |
 *  -------------------------------------------------
 */
    static constexpr real mu0 = 1.2566370614359172E-06;
    static constexpr real PI = 3.141592653589793;
    static const int NVTR_DIM = 5;  //!< Размерность массива nvtr.
    static const int STUFF_IND = 4; //!< Индекс номера материала в nvtr.
    static const int FE_DIM = 4;    //!< Кол-во узлов в конечном элементе.

    FileMap files;          //!< Карта файлов.
    int k;                  //!< Текущий конечный элемент.
    ActionType action;      //!< Текущее действие.
    real hx, hy;            //!< Текущие шаги по X и по Y.
    rvector2d G1;           //!< Первая матрица жесткости.
    rvector2d G2;           //!< Вторая матрица жесткости.
    rvector2d G3;           //!< Третья матрица жесткости.
    real eps;               //!< Точность решения.
    rvector solution;       //!< КЭ-решение задачи.
    intRealVector sources;  //!< Массив точечных источников.
    Slau matrix;

    FemPointPairList xySources;
    rvector powSources;
    FemPointPairList xyReceivers;
    intVector srcReceivers;
    rvector wReceivers;
    rvector errReceivers;

public: /*! Интерфейс. */
    FemSolver();
    FemSolver(FemDatGen &generator);
    FemSolver(FemDatGen &generator, const std::string & pointsFile, const std::string &path = ".");
    void input(const std::string & path = ".");
    void input(FemDatGen &generator, const std::string & path = ".");
    void inputPoints(const string &fname);
    void output();
    void solve();
    void solve(const FemStuffList &param);
    real value(real x, real y = 0);
    void setSigma(uint num, real value) {stuff.at(num).set_sig(value);}
    real sigma(uint num) {return stuff.at(num).sig();}
    void setPower(uint num, real value = 0) {powSources.at(num) = value;}
    //void setPower(real value) {powSources.front() = value;}
    real power(uint num = 0) {return powSources.at(num);}
    void setSourcePos(const FemPointPair & position){xySources.front() = position;}
    void setSourcePos(int num, const FemPointPair & position)
        {xySources.at(num) = position;}
    void test();
    int diffSolInPointsPair(rvector& sol);
    int diffSolInPointsPair0(rvector & sol);
    const rvector & getPointWeight() {return wReceivers;}
    const rvector & getPointNoise() {return errReceivers;}
    void getRadius(rvector & res);

private: /*! Закрытые функции. */
    void fillFileMap();
    void set_file(FileType ftype, const std::string & fname);
    std::string correct_path(const std::string & path);
    const std::string file(FileType ftype);

    void fillG();

    void getFromFile();
    void getNonlinearFiles();
    real _x(int i) {return xy[_global(i)].x;}
    real _y(int i) {return xy[_global(i)].y;}
    int _global(int i);
    void fillMatrix(Slau & matrix);
    void calc_h();

    real sig() {return stuff[nvtr[k][STUFF_IND]].sig();}
    bool point_lq_elem(real x, real y);

};
typedef FemSolver* FemSolver_ptr;

#endif // FEMSOLVER_H
