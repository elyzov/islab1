#ifndef FEMDATGEN_H
#define FEMDATGEN_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include "femnamespace.h"
#include "typedef.h"
#include "femarea.h"
#include "femgrid.h"
#include "femstuff.h"
#include "fempoint.h"
#include "femerror.h"

using namespace FEM;
/*! class FemDatGen - генерация .dat файлов для решателя МКЭ.
 *      ---------------------
 *      Формат входного файла:
 *      ---------------------
 *  KolObl
 *  x0[1]   x1[1]   y0[1]   y1[1]   nmat[1]
 *  x0[2]   x1[2]   y0[2]   y1[2]   nmat[2]
 *  ...
 *  x0[KolObl] x1[KolObl] y0[KolObl] y1[KolObl] nmat[KolObl]
 *
 *  KolMat
 *  nmat[1]   muo[1]   j[1]
 *  nmat[2]   muo[2]   j[2]
 *  ...
 *  nmat[KolMat]   muo[KolMat]   j[KolMat]
 *
 *  x0 KolX
 *  x[1]   hXmin[1]   dhX[1]   shX[1]
 *  x[2]   hXmin[2]   dhX[2]   shX[2]
 *  ...
 *  x[KolX] hXmin[KolX] dhX[KolX] shX[KolX]
 *
 *  y0 KolY
 *  y[1]   hYmin[1]   dhY[1]   shY[1]
 *  y[2]   hYmin[2]   dhY[2]   shY[2]
 *  ...
 *  y[KolY] hYmin[KolY] dhY[KolY] shY[KolY]
 *
 *  DoubleToX DoubleToY
 *
 *      ----------------------
 *      Формат выходных файлов:
 *      ----------------------
 *  inf2tr.dat (текстовый файл):
 *   kuzlov - число узлов сетки
 *   ktr    - число конечных элементов
 *   kt1    - число узлов с первыми краевыми условиями
 *
 *  nvtr.dat (двоичный файл):
 *   число записей - ktr
 *   структура i-й записи: 6*long (i1,i2,i3,i4,0,1),
 *   где i1,i2,i3,i4 - номера (глобальные) четырех вершин i-го прямоугольника
 *
 *  nvkat2d.dat (двоичный файл):
 *   число записей - ktr
 *   структура i-й записи: 1*long
 *   (номер материала i-го прямоугольника в соответствии с табл.1)
 *
 *  rz.dat (двоичный файл):
 *   число записей - kuzlov
 *   структура i-й записи: 2*double (x,y),
 *   где x,y - (x,y)-координаты i-й вершины
 *
 *  l1.dat (двоичный файл):
 *   число записей - kt1
 *   структура i-й записи: 1*long
 *   (номер i-й вершины с первым нулевыми краевым условием)
 *
 *  mu (текстовый файл)
 *   число строк - число материалов
 *   структура i-й строки: номер_материала  соответствующее_значение_мю
 *
 *  toku (текстовый файл)
 *   число строк - число материалов
 *   структура i-й строки: номер_материала  соответствующее_значение_тока
 */
class FemDatGen
{
private:
/*! -------------------------------------------------
 *                  Входные данные                  |
 *  -------------------------------------------------
 */

/*! Число прямоугольных областей,
 *  на которые разбита вся расчетная область. */
    int kolObl;

/*! Описание каждой из областей (всего kolObl). */
    FemAreaList areas;

/*! Число точечных источников. */
    int nSources;
    pointRealVector pointSrc;

/*! Описание известных материалов. */
    int kolMat;         //!< количество материалов.
    FemStuffList stuff; //!< описание материалов.

/*! Описание сетки по OX. */
    real x0;            //!< левая граница.
    int kolX;           //!< количество границ.
    FemGridList gridX;  //!< описание каждой границы (всего kolX).

/*! Описание сетки по OY. */
    real y0;            //!< левая граница.
    int kolY;           //!< количество границ.
    FemGridList gridY;  //!< описание каждой границы (всего kolY).

/*! Описание режима удвоения сетки по OX (DoubleToX) и OY (DoubleToY),
 *   0 - удвоения нет.
 *   1 - удвоение.
 *   2 - учетверение (двойное удвоение). */
    int doubleToX, doubleToY;


/*! -------------------------------------------------
 *                  Выходные данные                 |
 *  -------------------------------------------------
 */
    int kuzlov;         //!< число узлов сетки.
    int ktr;            //!< число конечных элементов.
    int kt1;            //!< число узлов с первыми краевыми условиями.

    intVector2d nvtr;   //!< массив конечных элементов (с номером материала).
    FemPointList xy;    //!< массив координат.
    intVector nvtrFst;  //!< массив первых краевых.
    intRealVector sources;   //!< массив точечных источников.

/*! -------------------------------------------------
 *              Вспомагательные данные              |
 *  -------------------------------------------------
 */
    FileMap files;          //!< Карта файлов.
    rvector vx, vy;         //!< points coordinates in OX, OY axis.
    int nx, ny;             //!< Кол-во точек по X и Y.
    int ndx, ndy;           //!< Кол-во отрезков по X и Y.

    static constexpr real cmp_eps = 1.e-4;  //!< Погрешность сравнения.
    static const int NVTR_DIM = 5;      //!< Размерность массива nvtr.
    static const int STUFF_IND = 4;     //!< Индекс номера материала в nvtr.

public: /*! интерфейс. */
    FemDatGen();
    FemDatGen(const std::string &file_in, const std::string &path_out, real zoom);
    FemDatGen(const std::string & file_in, real zoom = 0);
/*! Генерация. */
    void generate(const std::string & file_in,  //!< Файл с заданной областью.
                  const std::string & path_out, //!< Папка для вывода файлов.
                  real zoom = 0                 //!< Коэффицент масштабирования.
            );
    void generate(const std::string & file_in,  //!< Файл с заданной областью.
                  real zoom = 0                 //!< Коэффицент масштабирования.
            );
    void toTelmaGenerate(
            const std::string & fileIn,  //!< Файл с заданной областью.
            const std::string & fileOut, //!< Файл для Telma.
            real zoom = 0                //!< Коэффицент масштабирования.
            );
    void dataOutput(const std::string & path_out);
/*! Очистка. */
    void clear();

/*! Интерфейс вывода данных. */
    const intVector2d & finitElements() {return nvtr;}
    const FemPointList & coords() {return xy;}
    const intVector & firstConditions() {return nvtrFst;}
    const FemStuffList & materials() {return stuff;}
    const intRealVector & pointSources() {return sources;}

private: /*! закрытые функции. */
    void fillFileMap();
    void set_file(FileType ftype, const std::string & fname);
    std::string correct_path(const std::string & path);
    const std::string file(FileType ftype);

    void input(const std::string &fname, real koef = 0);

    int divAxis(real startCoord, const FemGridList & grid, rvector &v);
    void startGen();

    void output();
    void doubleTo(FemGridList &grid, int times);

    size_t findArea(const FemPoint &p);
    bool pointInArea(const FemPoint & p, int nArea);

    int getPointInd(const FemPoint & p);
};

#endif // FEMDATGEN_H
