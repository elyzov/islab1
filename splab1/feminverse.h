#ifndef FEMINVERSE_H
#define FEMINVERSE_H

#include "femsolver.h"
#include "typedef.h"

class FemInverse
{
    FemSolver_ptr dirTask;

    rvector expRes;
    uint M, N;
    uint maxIter;

    rvector u0;     //!< Вектор фиксированных значений.
    real alpha;     //!< Параметр Регуляризации.
    rvector noise;

    static constexpr real RELAX_EPS = 1.E-10;
    static constexpr real SOLVE_EPS = 1.E-8;
    static constexpr real GAMMA = 1.E-3;

public:
    FemInverse(FemSolver_ptr task);
    void solveLinear(const std::string & file);
    void solveNonLinear(const std::string & file);

private:
    inline real genErr();
    void getResp(const intVector & nParam, const rvector & param, rvector & resp);
    void getDiffResp(const intVector & nParam, const rvector & param, rvector & dRes);
    void getDiffResp(const rvector2d & err, const rvector & u, rvector & res);
    real calcJ(const rvector & dRes, const rvector &w);
    real calcJa(const rvector & dRes, const rvector &w, const rvector &u);
    real calcJ(const rvector2d & err, const rvector &u);
    real calcJa(const rvector2d & err, const rvector &u);
    real calcAlpha(real J, const rvector & u);
    void latexNonlinearOutput(const rvector & u_exp, const rvector & u, real nev,
                     real J, uint iter, real time, bool solved = true);
    void latexLinearOutput(const rvector & u, real nev);
};

#endif // FEMINVERSE_H
