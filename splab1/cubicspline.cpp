#include "cubicspline.h"
#include <math.h>
#include <algorithm>

CubicSpline::CubicSpline()
{
    n = 0;
}

void CubicSpline::sol_di_slau()
{
    for (int i = 1; i <= n; i++)
    {
        di[1][i] -= di[2][i-1] * di[0][i] / di[1][i-1];
        spline[i].c -= spline[i-1].c * di[0][i] / di[1][i-1];
    }
    spline[n].c /= di[1][n];
    for (int i = n-1; i >= 0; i--)
    {
        spline[i].c -= di[2][i] * spline[i+1].c;
        spline[i].c /= di[1][i];
    }
}

void CubicSpline::fill_matr()
{
    di[0][0] = 0;
    di[1][0] = 1;
    di[2][0] = 0;
    spline[0].c = 0;
    di[0][n] = 0;
    di[1][n] = 1;
    di[2][n] = 0;
    spline[n].c = 0;

    real h, h_nxt;

    h_nxt = spline[1].x - spline[0].x;

    for( int i = 1; i < n; i++ )
    {
        h = h_nxt;
        h_nxt = spline[i+1].x - spline[i].x;
        di[0][i] = h;
        di[1][i] = 2 * ( h + h_nxt );
        di[2][i] = h_nxt;
        spline[i].c = 3 * ( (spline[i+1].f - spline[i].f) / h_nxt -
                            (spline[i].f - spline[i-1].f) / h );
    }
}

void CubicSpline::calc_kof()
{
    for (int i = 0; i <= n; i++)
        spline[i].a = spline[i].f;

    spline[0].b = spline[0].d = 0;
    for( int i = 1; i <= n; i++ )
    {
        real h = spline[i].x - spline[i-1].x;
        spline[i].d = (spline[i].c - spline[i-1].c) / h / 3;
        spline[i].b = h*spline[i].c - h*h *spline[i].d +
                      (spline[i].f - spline[i-1].f) / h;
    }
}

int CubicSpline::build(const std::string &file)
{
    std::ifstream fin(file.c_str());
    if (!fin) return -1;
    //printf("Build spline from file %s\n", file.c_str());

    fin >> n;
    f.resize(n);
    spline.resize(n);
    for (int i = 0; i < n; i++)
    {
        real _x, _f;
        fin >> _x >> _f;
        if (fin.eof()) return -1;
        spline[i].f = _f;
        spline[i].x = _x;
    }

    std::sort(spline.begin(), spline.end());

    di.resize(3, rvector(n));
    n--;

    fill_matr();
    sol_di_slau();
    di.clear();
    calc_kof();

    return 0;
}

void CubicSpline::clear()
{
    f.clear();
    spline.clear();
    n = 0;
}

real CubicSpline::value(real x, ubyte derivative)
{
    if (n == 0) throw CubicSpline_Build_Error("n == 0");
                //return std::numeric_limits<real>::quiet_NaN();

    std::vector<tuple>::iterator curr;

    if (x <= spline.front().x)
        throw CubicSpline_Build_Error("x out of range (left)");
            //curr = spline.begin();
    else if (x >= spline.back().x)
        throw CubicSpline_Build_Error("x out of range (right)");
        // curr = spline.end()-1;
    else
    {
        size_t i = 0;
        size_t j = n;
        while (i + 1 < j)
        {
            size_t k = i + (j - i) / 2;
            if (x <= spline[k].x) j = k;
            else i = k;
        }
        curr = spline.begin() + j;
    }

    real dx = x - curr->x;
    switch (derivative)
    {
    case 0 : return curr->a + (curr->b + (curr->c + curr->d * dx) * dx) * dx;
    case 1 : return curr->b + (2 * curr->c + 3 * curr->d * dx) * dx;
    case 2 : return 2 * curr->c + 6 * curr->d * dx;
    case 3 : return 6 * curr->d;
    default : return 0;
    }
}

real CubicSpline::left_boundary()
{
    if (n == 0) return std::numeric_limits<real>::quiet_NaN();
    return spline.front().x;
}

real CubicSpline::right_boundary()
{
    if (n == 0) return std::numeric_limits<real>::quiet_NaN();
    return spline.back().x;
}
