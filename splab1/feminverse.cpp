#include "feminverse.h"
#include <math.h>
#include "denseslae.h"
#include "au.h"
#include <chrono>
#include "randomvariable.h"
#include <stdio.h>
#include <time.h>
#include "stopwatch.h"
#include "femdebug.h"

FemInverse::FemInverse(FemSolver_ptr task) : dirTask(task), maxIter(100)
{
    /*! Моделирование данных. */
    dirTask->solve();
    dirTask->diffSolInPointsPair(expRes);
    N = expRes.size();
    noise = dirTask->getPointNoise();

    //cout << "e = " << expRes << endl;
    RandomVariable::inst().setNormalError(noise, expRes);
    //cout << "e* = " << expRes << endl;

//    real expVar = RandomVariable::inst().sampleVariance(expRes);
//    for (uint i = 0; i < N; ++i) noise[i] = noise[i] ? noise[i] * expVar : 1;
}

void FemInverse::solveLinear(const string &file)
{
    if (file.empty()) return;
    std::ifstream fin(file.c_str());
    if (!fin) throw Fem_OpenFile_Error(file);
    fin.exceptions(std::ifstream::failbit |
                   std::ifstream::badbit  |
                   std::ifstream::eofbit);

    fin >> M;
    u0.resize(M);
    DenseSLAE slae(M);
    rvector u(M);
    rvector2d err(M, rvector(N));
    real curJ = 0;
    bool use_alpha = false;

    //const rvector & noise = dirTask->getPointNoise();

    fin >> use_alpha;
    for (uint i = 0; i < M; ++i)
    {
        fin >> u.at(i) >> u0.at(i);
    }

    Stopwatch::inst().start("linear");
    for (uint i = 0; i < M; ++i)
    {
        dirTask->diffSolInPointsPair0(err.at(i));
        //RandomVariable::inst().setNormalError(noise, err.at(i));
    }

    curJ = calcJ(err, u);
    alpha = use_alpha ? calcAlpha(curJ, u) : 0;

    for (uint i = 0; i < M; ++i)
    {
        for (uint j = 0; j < M; ++j)
        {
            slae.A(i, j) = i == j ? alpha : 0;
            for (uint k = 0; k < N; ++k)
                slae.A(i, j) += err[i][k] * err[j][k];
        }
        slae.F(i) = alpha * u0[i];
        for (uint k = 0; k < N; ++k)
            slae.F(i) += err[i][k] * expRes[k];
    }

    slae.solve();
    for (uint i = 0; i < M; ++i) u[i] = slae.X(i);

    curJ = calcJa(err, u);

    real execTime = Stopwatch::inst().stop("linear");
    latexLinearOutput(u, curJ);
    //cout << "Result u = " << u << endl;

}

void FemInverse::solveNonLinear(const string &file)
{
    if (file.empty()) return;
    std::ifstream fin(file.c_str());
    if (!fin) throw Fem_OpenFile_Error(file);
    fin.exceptions(std::ifstream::failbit |
                   std::ifstream::badbit  |
                   std::ifstream::eofbit);

    fin >> M;

    u0.resize(M);
    DenseSLAE slae(M);
    rvector u(M), err(N), u_next(M), du(M), u_exp(M);
    rvector2d derr(M, rvector(N));
    intVector nParam(M);
    real prevJ = 0, curJ = 0, betta, nev = 0;
    bool solved = false;
    bool invalidParam = false;
    bool use_alpha = false;
    rvector w(N);// = dirTask->getPointWeight();
    uint iterNum;

    fin >> use_alpha;
    for (uint i = 0; i < M; ++i)
    {
        fin >> nParam[i]>> u[i] >> u0[i];
    }

    /*! Извлекаем истинные значения параметров. */
    for (uint i = 0; i < M; ++i) u_exp[i] = dirTask->sigma(nParam[i]);

    /*! Вычисляем значения весов для приемников. */
    //getResp(nParam, u, noise, err);
    dirTask->getRadius(err);
    for (uint i = 0; i < N; ++i) w[i] = err[i];

    Stopwatch::inst().start("nonlinear");
    getDiffResp(nParam, u, err);
    //prevJ = calcJ(err, w);
    //alpha = use_alpha ? calcAlpha(prevJ, u) : 0;
    prevJ = calcJ(err, w);

    for (iterNum = 0; iterNum != maxIter; ++iterNum)
    {
        cerr << "InvIter #" << iterNum << ", J = " << prevJ
             << ", nev = " << nev << endl << flush;
        for (uint i = 0; i < M; ++i)
        {
            for (uint j = 0; j < M; ++j)
                du.at(j) = i == j ? 1.05 * u[j] : u[j];
            getDiffResp(nParam, du, derr.at(i));
        }

        /*! Заполняем СЛАУ. */
        alpha = 0;
        while (true)
        {
            for (uint i = 0; i < M; ++i)
            {
                for (uint j = 0; j < M; ++j)
                {
                    slae.A(i, j) = i == j ? alpha : 0;
                    for (uint k = 0; k < N; ++k)
                        slae.A(i, j) += pow(w[k], 2) * (derr[j][k] - err[k]) *
                                        (derr[i][k] - err[k]) / 0.025 / u[i] / u[j];
                }
                slae.F(i) = - alpha * (u[i] - u0[i]);
                for (uint k = 0; k < N; ++k)
                    slae.F(i) -= err[k] * pow(w[k], 2) *
                                 (derr[i][k] - err[k]) / 0.05 / u[i];
            }
            slae.solve();
            for (uint i = 0; i < M; i++)
            {
                u_next[i] = u[i] + slae.X(i);
                if (u_next[i] < 0 || u_next[i] > 10.) invalidParam = true;
            }
            cerr << "Test alpha = " << alpha << ", u = " << u_next << endl << flush;
            if (invalidParam)
            {
                //cerr << "There is invalid param in alpha" << endl << flush;
                invalidParam = false;
                alpha = alpha ? alpha * 5. : 10.e-8;
            }
            else break;
        }

        betta = 1;
        while (betta > RELAX_EPS)
        {
            for (uint i = 0; i < M; i++)
            {
                u_next[i] = u[i] + betta * slae.X(i);
                if (u_next[i] < 0) invalidParam = true;
            }
            //cerr << "Test betta u = " << u_next << endl << flush;
            if (invalidParam)
            {
                throw Fem_Solver_Error("There is invalid param in betta");
//                cerr << "There is invalid param in betta" << endl << flush;
//                invalidParam = false;
//                betta /= 2;
//                continue;
            }
            getDiffResp(nParam, u_next, err);
            curJ = calcJ(err, w);
            //cerr << "nev = " << nev << ", J = " << curJ << endl;
            if (fabs(curJ - prevJ) / prevJ < SOLVE_EPS)
            {
                solved = true;
                break;
            }
            else if (curJ > prevJ)
            {
                betta /= 2;
                continue;
            }
            else break;
        }
        if (betta < RELAX_EPS || solved) break;
        prevJ = curJ;
        for (uint i = 0; i < M; ++i) u[i] = u_next[i];
        rvector diffu(M);
        for (uint i = 0; i < M; ++i) diffu[i] = (u_next[i] - u_exp[i]);
        nev = norm(diffu) / norm(u_exp);
        cerr << "Current u = " << u << endl << flush;
    }
    real execTime = Stopwatch::inst().stop("nonlinear");
    latexNonlinearOutput(u_exp, u, nev, curJ, iterNum, execTime, solved);
    //cout << "Result u = " << u << endl;
}

void FemInverse::getResp(const intVector &nParam, const rvector &param, rvector &resp)
{
//    uint K = 10;
//    rvector tmp_resp(N);
//    resp.resize(N, 0);
    for (uint i = 0; i < M; ++i)
    {
        dirTask->setSigma(nParam.at(i), param.at(i));
    }
//    for (uint k = 0; k < K; ++k)
//    {
        dirTask->solve();
        dirTask->diffSolInPointsPair(resp);
//        for (uint i = 0; i < N; ++i) resp[i] += tmp_resp[i];
//    }
    //for (uint i = 0; i < N; ++i) resp[i] /= K;
}

void FemInverse::getDiffResp(const intVector &nParam, const rvector &param, rvector &dRes)
{
    rvector resp(N);
    getResp(nParam, param, resp);
    //RandomVariable::inst().setNormalError(noise, resp);
    for (uint i = 0; i < N; ++i) dRes[i] = (expRes[i] -resp[i]);
    //real svar = RandomVariable::inst().sampleVariance(dRes);
    //for (uint i = 0; i < N; ++i) dRes[i] /= svar;
}

real FemInverse::calcJ(const rvector &dRes, const rvector &w)
{
    real res = 0;
    for (uint i = 0; i < N; ++i) res += pow(w.at(i) * dRes.at(i), 2);
    return res;
}

real FemInverse::calcJa(const rvector &dRes, const rvector &w, const rvector & u)
{
    if (alpha)
    {
        real diff = 0;
        for (uint i = 0; i < M; ++i) diff += pow(u[i] - u0[i],2);
        return calcJ(dRes, w) + alpha * diff;
    }
    return calcJ(dRes, w);
}

void FemInverse::getDiffResp(const rvector2d & err, const rvector & u, rvector & res)
{
    for (uint i = 0; i < N; ++i)
    {
        res[i] = 0;
        for (uint j = 0; j < M; ++j)
        {
            res[i] += err[j][i] * u[j];
        }
        res[i] -= expRes[i];
    }
}

real FemInverse::calcJ(const rvector2d & err, const rvector & u)
{
    real res = 0;
    for (uint i = 0; i < N; ++i)
    {
        real sum = 0;
        for (uint j = 0; j < M; ++j) sum += err[j][i] * u[j];
        sum -= expRes[i];
        res += sum * sum;
    }
    return res;
}

real FemInverse::calcJa(const rvector2d &err, const rvector &u)
{
    if (alpha)
    {
        real diff = 0;
        for (uint i = 0; i < M; ++i) diff += pow(u[i] - u0[i],2);
        return calcJ(err, u) + alpha * diff;
    }
    return calcJ(err, u);
}

real FemInverse::calcAlpha(real J, const rvector &u)
{
    real diff = 0;
    for (uint i = 0; i < M; ++i) diff += pow(u[i] - u0[i],2);
    return J * GAMMA / diff;
}

void FemInverse::latexNonlinearOutput(const rvector &u_exp, const rvector &u, real nev, real J, uint iter, real time, bool solved)
{
    uint uSize = u.size();
    char str[256];
    char mark = solved ? ' ' : '*';
    if (uSize > 1)
    {
        sprintf(str, "\\multirow{%d}*{} & %e & %e%c & \\multirow{%d}{%e} & \\multirow{%d}{%e} & \\multirow{%d}{%d} & \\multirow{%d}{%f}\\\\\n" ,
                uSize, u_exp[0], u[0], mark, uSize, nev, uSize, J, uSize, iter, uSize, time);
        cout << str;
        for (uint i = 1; i < uSize; ++i)
        {
            sprintf(str, "                 & %e & %e%c &  &  &  &  \\\\\n", u_exp[i], u[i], mark);
            cout << str;
        }
    }
    else
    {
        sprintf(str, " & %e & %e%c & %e & %e & %d & %f\\\\\n" ,
                u_exp[0], u[0], mark, nev, J, iter, time);
        cout << str;
    }
}

void FemInverse::latexLinearOutput(const rvector &u, real nev)
{
    uint uSize = u.size();
    char str[256];
    if (uSize > 1)
    {
        sprintf(str, "\\multirow{%d}*{} & %e & \\multirow{%d}{%e}\\\\\n" ,
                uSize, u[0], uSize, nev);
        cout << str;
        for (uint i = 1; i < uSize; ++i)
        {
            sprintf(str, "                 & %e &  \\\\\n", u[i]);
            cout << str;
        }
    }
    else
    {
        sprintf(str, " & %e & %e\\\\\n" ,
                u[0], nev);
        cout << str;
    }
}
