#include "slau.h"

Slau::Slau(uint dim, real eps, int maxIter) : n(dim), e(eps), maxiter(maxIter)
{
}

Slau::~Slau()
{
}

void Slau::solve(SolveMethod m, rvector &sol)
{
    time_t start, end;
    time( &start ); //!< Засекаем время.
    switch(m)
    {
    case LOS_WF : los_without_fac(); break;
    case MSGN_WF : msgn_without_fac(); break;
    case MSGN_LU : lu_fac(); msgn(); break;
    case LOS_LU : lu_fac(); los(); break;
    case LOS_DI : diag_fac(); los(); break;
    case BCG_LU : lu_fac(); bcg(); break;
    case BCG_DI : diag_fac(); bcg(); break;
    case GMRES_1_LU : lu_fac(); gmres(1); break;
    case GMRES_2_LU : lu_fac(); gmres(2); break;
    case GMRES_3_LU : lu_fac(); gmres(3); break;
    case GMRES_4_LU : lu_fac(); gmres(4); break;
    case GMRES_5_LU : lu_fac(); gmres(5); break;
    case GMRES_6_LU : lu_fac(); gmres(6); break;
    case GMRES_10_LU : lu_fac(); gmres(10); break;
    }
    time( &end ); //!< Снимаем показания секундомера.
    stime = difftime(end, start);
    sol = x;
}

void Slau::mul_on_matrix(const rvector &v, rvector &res)
{
    mul_matrix(al, au, v, res);
}

real Slau::test_residual(const rvector &v)
{
    rvector res(n);
    mul_matrix(al, au, v, res);
//    out_vec(cout, "mul_res: ", res, n);
//    cout << endl << endl;
//    out_vec(cout, "f: ", f, n);
    for (int i = 0; i < n; i++) res[i] -= f[i];
    return sqrt(scalar(res, res) / scalar(f, f));
}

void Slau::mul_matrix(crvector &_al, crvector &_au, crvector &_x, rvector &res)
{
    for (uint i = 0; i < n; ++i)
    {
        res[i] = di[i] * _x[i];
        for (uint j = ia[i], end = ia[i+1]; j < end; ++j)
        {
            res[ ja[j] ] += _au[j] * _x[i];
            res[i] += _al[j] * _x[ ja[j] ];
        }
    }
}

void Slau::mul_up(crvector &_au, rvector &_x)
{
    if (_au.empty())
    {
        for (uint i = 0; i < n; ++i) _x[i] = gdi[i] * _x[i];
    }
    else
    {
        for (uint i = 0; i < n; ++i)
        {
            _x[i] = gdi[i] * _x[i];
            for (uint j = ia[i], end = ia[i+1]; j < end; ++j)
            {
                _x[ ja[j] ] += _au[j] * _x[i];
            }
        }
    }
}

void Slau::sol_up(crvector &_au, rvector &_f)
{
    if (_au.empty())
    {
        for (int i = n-1; i >= 0; --i) _f[i] /= gdi[i];
    }
    else
    {
        for (int i = n-1; i >= 0; --i)
        {
            _f[i] /= gdi[i];
            for (uint j = ia[i], end = ia[i+1]; j < end; ++j)
            {
                _f[ ja[j] ] -= _au[j] * _f[i];
            }

        }
    }
}

void Slau::sol_low( crvector &_al, rvector &_f )
{
    if (_al.empty())
    {
        for (uint i = 0; i < n; ++i) _f[i] /= gdi[i];
    }
    else
    {
        for (uint i = 0; i < n; ++i)
        {
            for (uint j = ia[i], end = ia[i+1]; j < end; ++j)
            {
                _f[i] -= _al[j] * _f[ ja[j] ];
            }
            _f[i] /= gdi[i];
        }
    }
}

void Slau::diag_fac()
{
    ggu.clear();
    ggl.clear();
    gdi.resize(n);
    for( int i = 0; i < n; i++ )
    {
        if( di[i] < 0 )
        {
            cerr << "Impossible to make the DIAG decomposition of the matrix." << endl;
            exit(1);
        }
        gdi[i] = sqrt( di[i] );
    }
}

void Slau::lu_fac()
{
    for (uint i = 0; i < n; ++i)
    {
        gdi[i] = di[i];
        for (uint j = ia[i], end = ia[i+1]; j < end; ++j)
        {
            real sum_l = 0;
            real sum_u = 0;
            uint k = ja[j];
            uint on_k = ia[k];
            uint on_i = ia[i];
            uint end_k = ia[k+1];
            while ( on_k < end_k && on_i < j )
            {
                if ( ja[on_i] == ja[on_k] )
                {
                    sum_l += ggl[on_i] * ggu[on_k];
                    sum_u += ggu[on_i] * ggl[on_k];
                    on_i++; on_k++;
                }
                else if ( ja[on_i] < ja[on_k] ) on_i++;
                else on_k++;
            }
            ggl[j] = ( al[j] - sum_l ) / gdi[k];
            ggu[j] = ( au[j] - sum_u ) / gdi[k];
            gdi[i] -= ggu[j] * ggl[j];
        }
        if ( gdi[i] < 0 )
        {
            cerr << "Impossible to make the LU(sq) decomposition of the matrix." << endl;
            exit(1);
        }
        gdi[i] = sqrt( gdi[i] );
    }
}

real Slau::scalar(crvector &v1, crvector &v2)
{
    real res = 0;
    for (uint i = 0; i < n; ++i) res += v1[i] * v2[i];
    return res;
}

void Slau::msgn_without_fac()
{
    rvector r(n);
    rvector z(n);
    rvector tarr(n);
    rvector tezz(n);
    real alpha, betta, rr, eps;

/*! Начальные значения. */
    it = 0;
//    out_vec(cout, "al: ", al, k);
//    out_vec(cout, "au: ", au, k);
//    out_vec(cout, "di: ", di, n);
//    out_vec(cout, "f: ", f, n);

    eps = scalar( f, f ) * e;
    mul_matrix( al, au, x, z );
    for( int i = 0; i < n; i++ ) z[i] = f[i] - z[i];
    mul_matrix( au, al, z, r );
    for( int i = 0; i < n; i++ ) z[i] = r[i];
    nev = scalar( r, r );

/*! Итерационный процесс. */
    while( nev > eps && it < maxiter )
    {
        it++;
        mul_matrix( al, au, z, tarr );
        mul_matrix( au, al, tarr, tezz );
        rr = nev;
        alpha =  rr / scalar( tezz, z );
        for( int i = 0; i < n; i++ )
        {
            x[i] += alpha * z[i];
            r[i] -= alpha * tezz[i];
        }
        nev = scalar( r, r );
        betta = nev / rr;
        for( int i = 0; i < n; i++ ) z[i] = r[i] + betta * z[i];
    }
}

void Slau::los_without_fac()
{
    rvector r(n);
    rvector z(n);
    rvector p(n);
    rvector tarr(n);
    real alpha, betta, pp;

/*! Начальные значения. */
    it = 0;

    mul_matrix( al, au, x, z );

    for( int i = 0; i < n; i++ )
    {
        r[i] = f[i] - z[i];
        z[i] = r[i];
    }
    //out_vec(cout, "z: ", z, n);

    mul_matrix( al, au, z, p );
    nev = scalar( r, r );

/*! Итерационный процесс. */
    while( nev > e && it < maxiter )
    {
        it++;
        pp =  scalar( p, p );
        if( pp == 0 )
        {
            cerr << "Div on 0 in iter = " << it << "." << endl;
            exit(1);
        }
        alpha = scalar( p, r ) / pp;
        for( int i = 0; i < n; i++ )
        {
            x[i] += alpha * z[i];
            r[i] -= alpha * p[i];
        }
        mul_matrix( al, au, r, tarr );
        betta = - scalar( p, tarr ) / pp;
        for( int i = 0; i < n; i++ )
        {
            z[i] = r[i] + betta * z[i];
            p[i] = tarr[i] + betta * p[i];
        }
        nev -= alpha * alpha * pp;
    }
}

void Slau::msgn()
{
    real alpha, betta, rr_next, rr, eps;

/*! Начальные значения. */
    it = 0;

    eps = sqrt(scalar( f, f )) * e;
    mul_matrix( al, au, x, z );
    for (uint i = 0; i < n; ++i)
    {
        r[i] = f[i] - z[i];
        tarr[i] = r[i];
    }
    sol_low(ggl, tarr);
    sol_up(ggu, tarr);
    nev = sqrt(scalar( r, r ));
    for (uint i = 0; i < n; ++i) z[i] = tarr[i];
    rr = scalar(tarr, r);

/*! Итерационный процесс. */
    while (nev > eps && it < maxiter)
    {
        ++it;
        mul_matrix( al, au, z, tezz );
        alpha = rr / scalar( tezz, z );
        for (uint i = 0; i < n; ++i)
        {
            x[i] += alpha * z[i];
            r[i] -= alpha * tezz[i];
            tarr[i] = r[i];
        }
        sol_low(ggl, tarr);
        sol_up(ggu, tarr);
        rr_next = scalar(tarr, r);
        nev = sqrt(scalar( r, r ));
        betta = rr_next / rr;
        for (uint i = 0; i < n; ++i) z[i] = tarr[i] + betta * z[i];
        std::swap(rr_next, rr);
    }
}

void Slau::los()
{
    rvector r(n);
    rvector z(n);
    rvector p(n);
    rvector tarr(n);
    rvector tezz(n);
    real alpha, betta, pp;

//    out_vec(cout, "al: ", al, k);
//    out_vec(cout, "au: ", au, k);
//    out_vec(cout, "di: ", di, n);
//    out_vec(cout, "f: ", f, n);
/*! Начальные значения. */
    it = 0;
    mul_matrix( al, au, x, z );
    for( int i = 0; i < n; i++ ) r[i] = f[i] - z[i];
    sol_low( ggl, r );
    for( int i = 0; i < n; i++ ) z[i] = r[i];
    sol_up( ggu, z );
    mul_matrix( al, au, z, p );
    sol_low( ggl, p );
    nev = scalar( r, r );

/*! Итерационный процесс. */
    while( nev > e && it < maxiter )
    {
        it++;
        pp = scalar( p, p );
        alpha = scalar( p, r ) / pp;
        for( int i = 0; i < n; i++ )
        {
            x[i] += alpha * z[i];
            r[i] -= alpha * p[i];
            tarr[i] = r[i];
        }
        sol_up( ggu, tarr );
        mul_matrix( al, au, tarr, tezz );
        sol_low( ggl, tezz );
        betta = - scalar( p, tezz ) / pp;
        for( int i = 0; i < n; i++ )
        {
            z[i] = tarr[i] + betta * z[i];
            p[i] = tezz[i] + betta * p[i];
        }
        nev -= alpha * alpha * pp;
    }
}

void Slau::bcg()
{
    rvector r(n);
    rvector z(n);
    rvector p(n);
    rvector s(n);
    rvector tarr(n);
    rvector tezz(n);
    real alpha, betta, pr, _pr;

/*! Начальные значения. */
    it = 0;
    mul_matrix(al, au, x, p);
    for (int i = 0; i < n; i++) r[i] = f[i] - p[i];
    sol_low(ggl, r);
    for (int i = 0; i < n; i++)
    {
        z[i] = p[i] = s[i] = r[i];
    }
    sol_up(ggu, z);
    nev = scalar( r, r );
    //real nev0 = nev;

    _pr = scalar(p, r);

/*! Итерационный процесс. */
    while( nev > e && it < maxiter )
    {
        it++;
        mul_matrix(al, au, z, tarr);
        sol_low(ggl, tarr);
        alpha = _pr / scalar(s, tarr);
        for (int i = 0; i < n; i++)
        {
            x[i] += alpha * z[i];
            r[i] -= alpha * tarr[i];
            tezz[i] = s[i];
        }
        sol_up(ggl, tezz);
        mul_matrix(au, al, tezz, tarr);
        sol_low(ggu, tarr);
        for (int i = 0; i < n; i++)
        {
            p[i] -= alpha * tarr[i];
            tezz[i] = r[i];
        }
        pr = scalar(p, r);
        betta = pr / _pr;
        sol_up(ggu, tezz);
        for (int i = 0; i < n; i++)
        {
            z[i] = tezz[i] + betta * z[i];
            s[i] = p[i] + betta * s[i];
        }
        swap(pr, _pr);
        nev = scalar(r,r);
    }
}

void Slau::gmres(int m)
{
    rvector r(n);
    rvector w(n);
    rvector tarr(n);
    rvector2d H(m+1, rvector(m));
    rvector d(m+1);
    rvector2d v(n, rvector(m));
//    alloc_arr(H, m+1, m);
//    alloc_arr(v, n, m);

/*! Начальные значения. */
    it = 0;
    mul_matrix(al, au, x, r);
    mul_up(ggu, x);
    for (int i = 0; i < n; i++) r[i] = f[i] - r[i];
    sol_low(ggl, r);
    nev = sqrt(scalar(r, r));

/*! Итерационный процесс. */
    while (nev > e && it < maxiter)
    {
        //cout << it << flush;
        for (int j = 0; j < n; j++) v[j][0] = r[j] / nev;

        for (int u = 0; u < m; u++)
        {
            for (int i = 0; i < n; i++) tarr[i] = v[i][u];
            sol_up(ggu, tarr);
            mul_matrix(al, au, tarr, w);
            sol_low(ggl, w);

            for (int l = 0; l <= u; l++)
            {
                H[l][u] = 0;
                for (int i = 0; i < n; i++) H[l][u] += v[i][l] * w[i];
                for (int i = 0; i < n; i++) w[i] -= H[l][u] * v[i][l];
            }
            for (int l = u+1; l <= m; l++) H[l][u] = 0;

            real _H = H[u+1][u] = sqrt(scalar(w, w));
            if (_H == 0)
            {
                m = u + 1;
                break;
            }
            if (u != m-1)
            {
                for (int i = 0; i < n; i++) v[i][u+1] = w[i] / _H;
            }
        }

        d[0] = nev;
        for (int i = 1; i <= m; i++) d[i] = 0;

//        H[0][0] = 1, H[0][1] = 2, H[0][2] = 3;
//        H[1][0] = 1, H[1][1] = 4, H[1][2] = 5;
//        H[2][0] = 0, H[2][1] = 1, H[2][2] = 6;

//        d[0] = 6, d[1] = 10, d[2] = 7;

//        for (int i =0; i < m; i++) out_vec(cout, "H: ", H[i], m);

//        for (int i = 1; i < m; i++)
//        {
//            real _H = H[i][i-1] / H[i-1][i-1];
//            for (int j = i; j < m; j++)
//                H[i][j] -= H[i-1][j] * _H;
//            d[i] -= d[i-1] * _H;
//        }

        for (int i = 0; i < m; i++)
        {
            real s = H[i+1][i] / sqrt(pow(H[i][i], 2) + pow(H[i+1][i], 2));
            real c = H[i][i] / sqrt(pow(H[i][i], 2) + pow(H[i+1][i], 2));
            for (int j = i; j < m; j++)
            {
                real _H_i = H[i][j];
                real _H_i1 = H[i+1][j];
                real _d_i = d[i];
                real _d_i1 = d[i+1];

                H[i][j] = c * _H_i + s * _H_i1;
                H[i+1][j] = -s * _H_i + c * _H_i1;
                d[i] = c * _d_i + s * _d_i1;
                d[i+1] = -s * _d_i + c * _d_i1;
            }
        }

        for (int i = m-1; i >= 0; i--)
        {
            for (int j = m-1; j > i; j--)
                d[i] -= d[j] * H[i][j];
            d[i] /= H[i][i];
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
                x[i] += d[j] * v[i][j];
        }

        for (int i = 0; i < n; i++) tarr[i] = x[i];
        sol_up(ggu, tarr);
        mul_matrix(al, au, tarr, r);
        for (int i = 0; i < n; i++) r[i] = f[i] - r[i];
        sol_low(ggl, r);
        nev = sqrt(scalar(r, r));
        it++;
    }
    sol_up(ggu, x);
}

void Slau::gen_with_nvtr(const intVector2d &_nvtr, int _n, int _m)
{
    ilvector l(n);
    intVector2d nvtr(_n, intVector(_m));

/*! Сотрировка элементов строк массива. */
    ilist for_sort;
    for( int i = 0; i < _n; i++ )
    {
        for( int j = 0; j < _m; j++ )
            for_sort.push_back(_nvtr[i][j]);
        for_sort.sort();

        for( int j = 0; j < _m; j++ )
        {
            nvtr[i][j] = for_sort.front();
            for_sort.pop_front();
        }
    }
    for_sort.clear();

//    for (int i = 0; i < _n; i++)
//        out_vec(cout, "nvtr 2: ", nvtr[i], _m);

/*! Заполняем массив списков. */
    for( int i = 0; i < _n; i++ )
        for( int j = 1; j < _m; j++ )
            for( int k = j - 1; k >= 0; k-- )
                if( nvtr[i][j] > nvtr[i][k] )
                    l[ nvtr[i][j] ].push_back( nvtr[i][k] );
    nvtr.clear();

/*! Сортируем и удаляем повторяющиеся элементы. */
    for( int i = 0; i < n; i++ ){ l[i].sort(); l[i].unique(); }

//    for( int i = 0; i < n; i++ )
//    {
//        cout << "list #" << i+1 << ": ";
//        for( _i_list j = l[i].begin(); j != l[i].end(); ++j )
//            cout << *j << ' ';
//        cout << endl;
//    }

/*! Подготавливаем массив ia. */
    ia.resize(n+1);
    ia[0] = 0;

/*! Заполняем массив ia. */
    for( int i = 0; i < n; i++ ) ia[ i+1 ] = ia[i] + l[i].size();

/*! Подготавливаем массив ja. */
    k = ia[n] - ia[0];
    ja.resize(k);

/*! Заполяем массив ja. */
    int _k = 0;
    for( int i = 0; i < n; i++ )
    {
        ia[ i+1 ] = ia[i] + l[i].size();
        for( _ilist it = l[i].begin(); it != l[i].end(); ++it )
        {
            ja[_k] = *it; _k++;
        }
    };

    au.resize(k);
    al.resize(k);
    ggu.resize(k);
    ggl.resize(k);
    di.resize(n);
    gdi.resize(n);
    x.resize(n);
    f.resize(n);

    r.resize(n);
    z.resize(n);
    tarr.resize(n);
    tezz.resize(n);
}

real & Slau::A(int i, int j)
{
    if( i == j ) return di[i];
    bool is_up;
    if (is_up = (i < j)) swap(i, j);

    for( int k = ia[i], end = ia[i+1]; k < end; k++ )
    {
        if( ja[k] == j ) return is_up ? au[k] : al[k];
    }
    cerr << "Не найден элемент (" << i << ',' << j << ") в матрице A." << endl;
    exit(1);
}

real Slau::__A(int i, int j)
{
    if( i == j ) return di[i];
    bool is_up;
    if (is_up = (i < j)) swap(i, j);

    for( int k = ia[i], end = ia[i+1]; k < end; k++ )
    {
        if( ja[k] == j ) return is_up ? au[k] : al[k];
    }
    return 0;
}

void Slau::set_first(int i, real val)
{
    di[i] = 1.;
    f[i] = val;
    for (uint k = ia[i], end = ia[i+1]; k < end; ++k)
    {
        f[ ja[k] ] -= au[k] * val;
        au[k] = 0.;
        al[k] = 0.;
    }
    for (uint on_i = i+1; on_i < n; ++on_i)
    {
        for (uint k = ia[on_i], end = ia[on_i+1]; k < end; ++k)
        {
            if( ja[k] == i )
            {
                f[ on_i ] -= al[k] * val;
                al[k] = 0.;
                au[k] = 0.;
            }
        }
    }
}

void Slau::init(uint dim, real eps, int maxIter)
{
    n = dim;
    e = eps;
    maxiter = maxIter;
}

void Slau::clear()
{
    for (uint i = 0; i < k; ++i)
        { al[i] = 0; au[i] = 0; }
    for (uint i = 0; i < n; ++i)
        { di[i] = 0; x[i] = 0; f[i] = 0; }
}

void Slau::print()
{
    out_vec(cout, "ia: ", ia, n+1);
    out_vec(cout, "ja: ", ja, k);
    out_vec(cout, "al: ", al, k);
    out_vec(cout, "au: ", au, k);
    out_vec(cout, "di: ", di, n);
    out_vec(cout, "f: ", f, n);
    out_vec(cout, "ggl: ", ggl, k);
    out_vec(cout, "ggu: ", ggu, k);
    out_vec(cout, "gdi: ", gdi, n);
}

void Slau::metrixForm()
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << __A(i,j) << " ";
        cout << endl;
    }
}
