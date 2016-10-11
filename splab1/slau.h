#ifndef SLAU_H
#define SLAU_H

#include <iostream>
#include <fstream>
#include <list>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "typedef.h"

using namespace std;

/*! class Slau - решение СЛАУ. Матрица хранится в разреженном
 *  строчно-столбцовом формате. Портрет матрицы симметричный. */
class Slau
{
    uint n, k;           //!< [in] Размерность матрицы, ко-во элементов.
    intVector ia, ja;   //!< [in] Индаксация по столбцам/строкам.
    rvector al, au;     //!< [in] Нижний/верхний треугольники.
    rvector di;         //!< [in] Диагональ.
    rvector f, x;       //!< [in] Правая часть, вектор решения.

/*! Матрица предобуславливания. */
    rvector ggu, ggl;   //!< Нижний/верхний треугольники.
    rvector gdi;        //!< Диагональ.

    rvector r, z, tarr, tezz;

    real e;             //!< [in] Точность вычислений.
    uint maxiter;        //!< [in] Максимальное кол-во итераций.

    real nev;           //!< [out] Невязка.
    uint it;             //!< [out] Кол-во итерации.
    int style;          //!< [in/own] Стиль задания массивов.
    double stime;       //!< [out] Время решения задачи.

public: /*! Интерфейс класса. */
    enum SolveMethod { LOS_WF, MSGN_WF, MSGN_LU, LOS_LU, LOS_DI,
                    BCG_LU, BCG_DI, GMRES_3_LU, GMRES_6_LU,
                    GMRES_10_LU, GMRES_1_LU, GMRES_2_LU,
                    GMRES_4_LU, GMRES_5_LU };

    Slau(uint dim = 0, real eps = 1.e-7, int maxIter = 20000);
    ~Slau();

/*! SET-методы. */
    void set_e(real eps){ e = eps; }
    void set_first( int i, real val );
    void dub_al(){ al = au; }
    void dub_au(){ au = al; }
    void init(uint dim, real eps = 1.e-7, int maxIter = 20000);

/*! GET-методы. */
    const rvector & solution() { return x; }
    const rvector & rightPart() { return f; }
    real residual(){ return nev; }
    double solutionTime(){ return stime; }
    int solutionIterations(){ return it; }
    real & A( int i, int j );       //!< i-ый, j-ый элемент матрицы.
    real & F(int i){ return f[i]; } //!< i-ый элемент правой части.

/*! Методы ввода. */
    void gen_with_nvtr(const intVector2d &nvtr, int _n, int _m);

/*! Методы вывода. */
    void metrixForm();

/*! Решение СЛАУ. */
    void solve(SolveMethod m, rvector & sol);

/*! Методы для отладки. */
    void mul_on_matrix(const rvector &v, rvector &res);
    real test_residual(const rvector &v);
    void print();

/*! Обнуление СЛАУ. */
    void clear();


private:
    real __A( int i, int j );   //!< i-ый, j-ый элемент матрицы.
/*! Умножение матрицы на вектор. */
    void mul_matrix(crvector &_al, crvector &_au, crvector &_x, rvector &res);

/*! Умножение верхне-треугольной матрицы на вектор. */
    void mul_up(crvector &_au, rvector &_x);

/*! Решение верхне-треугольной СЛАУ. */
    void sol_up(crvector &_au, rvector &_f);

/*! Решение нижне-треугольной СЛАУ. */
    void sol_low(crvector &_al, rvector &_f);

/*! Диагональная факторизация. */
    void diag_fac();

/*! LU(sq) факторизация. */
    void lu_fac();

/*! Вычисление скалярного произведения. */
    real scalar(crvector &v1, crvector &v2);

/*! МСГ без факторизации. */
    void msgn_without_fac();

    void msgn();

/*! ЛОС без факторизации. */
    void los_without_fac();

/*! ЛОС с факторизацией. */
    void los();

/*! BCG с факторизацией. */
    void bcg();

/*! GMRES. */
    void gmres(int m);
};

/*! Тестовый вывод. */
template <typename type>
void out_vec( ostream & cout, const char * name, const type & arr, int n )
{
    cout << name << ' ';
    if( !arr.empty() ) for( int i = 0; i < n; i++ ) cout << arr[i] << ' ';
    else cout << "is empty";
    cout << endl;
}

#endif // SLAU_H
