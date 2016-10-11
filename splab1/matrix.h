#ifndef MATRIX_H
#define MATRIX_H

#include "typedef.h"


template <class Type>
class Matrix
{
    uint n; //!< Row count.
    uint m; //!< Col count.
    std::vector<Type> data;
public:
    Matrix(uint rowCount, uint colCount) : n(rowCount), m(colCount)
    {
        data.resize(n * m);
    }

    Matrix(){}

    Type & operator ()(uint row, uint col)
    {
        return data.at(row * m + col);
    }

    void resize(uint rowCount, uint colCount)
    {
        n = rowCount;
        m = colCount;
        data.resize(n * m);
    }
    void swapColumns(uint i, uint j)
    {
        for (uint k = 0; k < n; ++k) std::swap(matr(k,i), matr(k,j));
    }

    uint dimH();
    uint dimV();

private:
    Type & matr(uint row, uint col)
    {
        return data.at(row * m + col);
    }

};

#endif // MATRIX_H
