#include "denseslae.h"

DenseSLAE::DenseSLAE(uint dim) : _dim(dim), a(dim, dim), f(dim), x(dim)
{
}

void DenseSLAE::solve()
{
    //if (_dim == 1) x[0] = f[0] / a(0, 0);
    uintVector mask(_dim);   //!< Маска востановления индексов неизвестных.
    for (uint i = 0; i < _dim; ++i) mask[i] = i;

    /*! Прямой ход. */
    for (uint i = 0; i < _dim; ++i)
    {
        uint maxInd = i;
        real maxVal = a(i,i);

        /*! Поиск главного элемента в строке. */
        for (uint j = i+1; j < _dim; ++j)
            if (fabs(maxVal) < fabs(a(i,j)))
            {
                maxVal = a(i,j);
                maxInd = j;
            }

        if (maxVal == 0)
            throw DenseSLAE_Error("Can't solve SLAE: there is linearly dependence");

        /*! Перестановка столбцов. */
        if (i != maxInd)
        {
            a.swapColumns(i, maxInd);
            std::swap(mask[i], mask[maxInd]);
        }

        /*! Сам ход. */
        for (uint j = i; j < _dim; ++j) a(i, j) /= maxVal;
        f[i] /= maxVal;

        for (uint j = i+1; j < _dim; ++j)
        {
            for (uint k = i+1; k < _dim; ++k)
                a(j,k) -= a(j,i) * a(i,k);
            f[j] -= a(j,i) * f[i];
        }
    }

    /*! Обратный ход. */
    for (int i = _dim - 1; i >= 0; --i)
        for (int j = _dim - 1; j > i; --j)
            f[i] -= a(i,j) * f[j];

    /*! Восстанавливаем индексы. */
    for (uint i = 0; i < _dim; ++i) x[mask[i]] = f[i];
}

void DenseSLAE::test()
{
    for (uint i = 0; i < _dim; ++i)
    {
        f[i] = 0;
        for (uint j = 0 ; j < _dim; ++j)
        {
            a(i, j) = i == j ? 2 : 1;
            f[i] += a(i, j) * j;
        }
    }
}
