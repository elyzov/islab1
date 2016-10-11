#include "femsolver.h"
#include <stdio.h>
#include <math.h>

FemSolver::FemSolver()
{
    fillFileMap();
    fillG();
    eps = 1.e-6;
}

FemSolver::FemSolver(FemDatGen &generator)
{
    fillFileMap();
    fillG();
    eps = 1.e-6;
    input(generator);
}

FemSolver::FemSolver(FemDatGen &generator, const std::string & pointsFile, const std::string &path)
{
    fillFileMap();
    fillG();
    eps = 1.e-6;
    input(generator, path);
    inputPoints(pointsFile);
}

void FemSolver::fillFileMap()
{
    files[FT_SIZE] = "inf2tr.dat";
    files[FT_NVTR] = "nvtr.dat";
    files[FT_NVTRSTUFF] = "nvkat2d.dat";
    files[FT_XY] = "rz.dat";
    files[FT_NVTRFST] = "l1.dat";
    files[FT_MU] = "mu";
    files[FT_J] = "toku";
    files[FT_RES] = "v2.dat";
    files[FT_WORKDIR] = ".";
    files[FT_B_RES] = "b_res.dat";
}

void FemSolver::fillG()
{
    G1.resize(FE_DIM, rvector(FE_DIM));
    G2.resize(FE_DIM, rvector(FE_DIM));
    G3.resize(FE_DIM, rvector(FE_DIM));
    for (int i = 0; i < FE_DIM; i++)
    {
        for (int j = 0; j < FE_DIM; j++)
        {
            if (i == j) {G1[i][j] = 2; G2[i][j] = 2;}
            else if (i + j == 3) {G1[i][j] = -1; G2[i][j] = -1;}
            else if (i + j == 1 || i + j == 5) {G1[i][j] = -2; G2[i][j] = 1;}
            else if (i + j == 2 || i + j == 4) {G1[i][j] = 1; G2[i][j] = -2;}
        }
    }
    G3[1][1] = G3[3][3] = 3.;
    G3[3][1] = G3[1][3] = -3.;
    G3[0][0] = G3[0][1] = G3[1][0] = G3[2][2] = G3[2][3] = G3[3][2] = 1.;
    G3[2][0] = G3[2][1] = G3[3][0] = G3[0][2] = G3[0][3] = G3[1][2] = -1.;

}

int FemSolver::_global(int i)
{
    switch (action)
    {
    case AT_LOCAL : return nvtr[k][i];
    case AT_FIRST : return i;
    default : throw Fem_Solver_Error("Unknown current action");
    }
}

void FemSolver::input(const std::string &path)
{
    set_file(FT_WORKDIR, path);
    getFromFile();
    getNonlinearFiles();
}

void FemSolver::input(FemDatGen &generator, const string &path)
{
    set_file(FT_WORKDIR, path);
    nvtr = generator.finitElements();
    xy = generator.coords();
    nvtrFst = generator.firstConditions();
    stuff = generator.materials();
    sources = generator.pointSources();
    kuzlov = xy.size();
    ktr = nvtr.size();
    kt1 = nvtrFst.size();
    getNonlinearFiles();
    matrix.init(kuzlov);
    matrix.gen_with_nvtr(nvtr, ktr, NVTR_DIM);
    matrix.set_e(1.e-15);
    solution.resize(kuzlov, 0);
}

void FemSolver::inputPoints(const string &fname)
{
    std::ifstream fin(fname.c_str());
    if (!fin) throw Fem_OpenFile_Error(fname);
    fin.exceptions(std::ifstream::failbit |
                   std::ifstream::badbit  |
                   std::ifstream::eofbit);

    uint num;
    fin >> num;
    xySources.resize(num);
    powSources.resize(num);
    for (uint i = 0; i < num; ++i)
    {
        fin >> xySources.at(i).first.x >> xySources.at(i).first.y;
        fin >> xySources.at(i).second.x >> xySources.at(i).second.y;
        fin >> powSources.at(i);
    }

    fin >> num;
    xyReceivers.resize(num);
    wReceivers.resize(num);
    srcReceivers.resize(num);
    errReceivers.resize(num);
    for (uint i = 0; i < num; ++i)
    {
        fin >> xyReceivers.at(i).first.x
            >> xyReceivers.at(i).first.y;
        fin >> xyReceivers.at(i).second.x
            >> xyReceivers.at(i).second.y;
        fin >> srcReceivers.at(i);
        fin >> wReceivers.at(i);
        fin >> errReceivers.at(i);
        errReceivers.at(i) *= 0.01;
    }
}

void FemSolver::output()
{
    std::ofstream fout;
    std::ios_base::openmode binary_mode = std::ios_base::out |
                                          std::ios_base::trunc |
                                          std::ios_base::binary;
    fout.open(file(FT_RES).c_str(), binary_mode);
    if (!fout) throw Fem_OpenFile_Error(file(FT_RES));
    for (int i = 0, end = kuzlov; i != end; i++)
    {
        double data = solution.at(i);
        fout.write((char *)&data, sizeof(double));
    }
    fout.close();
}

void FemSolver::solve(const FemStuffList &param)
{
    for (ciFemStuffList it = param.begin(); it != param.end(); ++it)
    {
        stuff.at(it->first) = it->second;
    }
    solve();
}

void FemSolver::solve()
{
    fillMatrix(matrix);
    matrix.solve(Slau::MSGN_LU, solution);
    //real nev = matrix.test_residual(solution);
}

void FemSolver::test()
{
    real C = 0.5 / PI;
    for (size_t i = 0; i < solution.size(); ++i)
    {
        real diff = fabs(solution.at(i) - C / xy.at(i).x) / C / xy.at(i).x;
        //if (diff > 1.E-4)
            cout << xy.at(i).x << ", " << xy.at(i).y << " : " << diff
                 << " : " << solution.at(i) << "* -> " << C / xy.at(i).x << endl;
    }

}

void FemSolver::fillMatrix(Slau &matrix)
{
    matrix.clear();
/*! Add local matrix's in global matrix. */
    action = AT_LOCAL;
    for (k = 0; k < ktr; k++)
    {
        calc_h();
        real _sig = sig();
        real _hc2 = _x(0) * hy / hx / 6.;
        real _hc3 = hy / 12.;
        real _hc1 = _hc2 + _hc3;
        for (uint i = 0; i < FE_DIM; ++i)
        {
            for (uint j = 0; j <= i; ++j)
            {
                real _G =  G1[i][j] * _hc1 + G2[i][j] * _hc2 + G3[i][j] * _hc3;
                matrix.A(_global(i),_global(j)) += _sig * _G;
            }
        }
    }
    matrix.dub_au();

/*! Set sources' power in right part. */
    for (size_t i = 0; i < sources.size(); i++)
    {
        matrix.F(sources.at(i).first) += sources.at(i).second;
    }

/*! Set first bounbary conditions. */
    action = AT_FIRST;
    for (iintVector it = nvtrFst.begin(); it != nvtrFst.end(); it++)
    {
        matrix.set_first(*it, 0.);
    }

/*! Second boundary conditions set automatically. */
}

void FemSolver::calc_h()
{
    for( int i = 1; i < FE_DIM; i++ )
        if( ( hx = fabs( _x(i) - _x(0) ) ) != 0 ) break;
    for( int i = 1; i < FE_DIM; i++ )
        if( ( hy = fabs( _y(i) - _y(0) ) ) != 0 ) break;
}

void FemSolver::getFromFile()
{
    std::ifstream fin;
    std::ios_base::openmode binary_mode = std::ios_base::in |
                                          std::ios_base::binary;
//    fin.exceptions(std::ifstream::failbit |
//                   std::ifstream::badbit  |
//                   std::ifstream::eofbit);

/*! Input sizes. */
    fin.open(file(FT_SIZE).c_str());
    if (!fin) throw Fem_OpenFile_Error(files[FT_SIZE]);
    fin >> kuzlov >> ktr >> kt1;
    fin.close();

    nvtr.resize(ktr, intVector(NVTR_DIM));
    xy.resize(kuzlov);
    nvtrFst.resize(kt1);

/*! Input nvtr. */
    fin.open(file(FT_NVTR).c_str(), binary_mode);
    if (!fin) throw Fem_OpenFile_Error(files[FT_NVTR]);
    for (int i = 0; i != ktr; i++)
    {
        long data[4];
        fin.read((char *)data, 4 * sizeof(long));
        nvtr[i][2] = (int)data[0] - 1;
        nvtr[i][3] = (int)data[1] - 1;
        nvtr[i][0] = (int)data[2] - 1;
        nvtr[i][1] = (int)data[3] - 1;
        fin.seekg(2 * sizeof(long), std::ios_base::cur);
    }
    fin.close();

/*! Input stuff for finite elements. */
    fin.open(file(FT_NVTRSTUFF).c_str(), binary_mode);
    if (!fin) throw Fem_OpenFile_Error(files[FT_NVTRSTUFF]);
    for (int i = 0; i != ktr; i++)
    {
        long data;
        fin.read((char *)&data, sizeof(long));
        nvtr[i][STUFF_IND] = (int)data;
    }
    fin.close();

/*! Input coords of points. */
    fin.open(file(FT_XY).c_str(), binary_mode);
    if (!fin) throw Fem_OpenFile_Error(files[FT_XY]);
    for (int i = 0; i != kuzlov; i++)
    {
        double data[2];
        fin.read((char *)data, 2*sizeof(double));
        xy[i].x = data[0]; xy[i].y = data[1];
    }
    fin.close();

/*! Input first boundary conditions. */
    fin.open(file(FT_NVTRFST).c_str(), binary_mode);
    if (!fin) throw Fem_OpenFile_Error(files[FT_NVTRFST]);
    for (int i = 0; i != kt1; i++)
    {
        long data;
        fin.read((char *)&data, sizeof(long));
        nvtrFst[i] = data - 1;
    }
    fin.close();

/*! Input mu0 from stuff. -- TEMPORARILY! WITHOUT ERROR PROCESS!*/
    fin.open(file(FT_MU).c_str());
    if (!fin) throw Fem_OpenFile_Error(files[FT_MU]);
    while (!fin.eof())
    {
        int nmat;
        real muo;
        fin >> nmat >> muo;
        stuff[nmat] = FemStuff(muo);
    }
    fin.close();

///*! Input J from stuff. -- TEMPORARILY! WITHOUT ERROR PROCESS!*/
//    fin.open(file(FT_J).c_str());
//    if (!fin) throw Fem_OpenFile_Error(files[FT_J]);
//    while (!fin.eof())
//    {
//        int nmat;
//        real j;
//        fin >> nmat >> j;
//        stuff[nmat].set_j(j);
//    }
//    fin.close();
}

void FemSolver::getNonlinearFiles()
{
//    char fname[7];
//    std::string toOpen;
//    for (iFemStuffList it = stuff.begin(); it != stuff.end(); ++it)
//    {
//        sprintf(fname, "mu.%03d", it->first);
//        toOpen = files[FT_WORKDIR] + fname;
//        it->second.test_nonlinearty(toOpen);
//    }
}

int FemSolver::diffSolInPointsPair(rvector & sol)
{
    if (sol.empty()) sol.resize(xyReceivers.size());
    for (uint i = 0, iend = xyReceivers.size(); i != iend; ++i)
    {
        FemPoint & pM = xyReceivers.at(i).first;
        FemPoint & pN = xyReceivers.at(i).second;
        FemPoint & pA = xySources.at(srcReceivers.at(i)).first;
        FemPoint & pB = xySources.at(srcReceivers.at(i)).second;
        sol.at(i) = value(len(pA, pN)) - value(len(pA, pM)) -
                    value(len(pB, pN)) + value(len(pB, pM));
        sol.at(i) *= powSources.at(srcReceivers.at(i));
    }
    return 0;
}

int FemSolver::diffSolInPointsPair0(rvector & sol)
{
    if (sol.empty()) sol.resize(xyReceivers.size());
    for (uint i = 0, iend = xyReceivers.size(); i != iend; ++i)
    {
        FemPoint & pM = xyReceivers.at(i).first;
        FemPoint & pN = xyReceivers.at(i).second;
        FemPoint & pA = xySources.at(srcReceivers.at(i)).first;
        FemPoint & pB = xySources.at(srcReceivers.at(i)).second;
        sol.at(i) = value(len(pA, pN)) - value(len(pA, pM)) -
                    value(len(pB, pN)) + value(len(pB, pM));
    }
    return 0;
}

void FemSolver::getRadius(rvector &res)
{
    if (res.empty()) res.resize(xyReceivers.size());
    for (uint i = 0, iend = xyReceivers.size(); i != iend; ++i)
    {
        FemPoint & pM = xyReceivers.at(i).first;
        FemPoint & pN = xyReceivers.at(i).second;
        FemPoint & pA = xySources.at(srcReceivers.at(i)).first;
        FemPoint & pB = xySources.at(srcReceivers.at(i)).second;
        res.at(i) = len(mid(pA, pB), mid(pN, pM));
    }
}

bool FemSolver::point_lq_elem(real x, real y)
{
    if (y < _y(0)) return true;
    if (y > _y(2)) return false;
    if (x < _x(0)) return true;
    if (x > _x(1)) return false;
    return true;
}

real FemSolver::value(real x, real y)
{
    if (solution.empty()) return std::numeric_limits<real>::quiet_NaN();

    if (x < xy.front().x || x > xy.back().x)
        throw Fem_Solver_Error("value out of range");
    if (y < xy.front().y || y > xy.back().y)
        throw Fem_Solver_Error("value out of range");

    action = AT_LOCAL;

    {   //!< Find finit element with point (X,Y)
        size_t i = 0;
        size_t j = ktr-1;
        while (i + 1 < j)
        {
            k = i + (j - i) / 2;
            if (point_lq_elem(x, y))
                j = k;
            else i = k;
        }
        k = j;
    }

//    __PRINT__(_x(0)); __PRINT__(_x(1));
//    __PRINT__(_y(0)); __PRINT__(_y(2));

    calc_h();
    return solution[_global(0)] * (_x(1) - x) * (_y(2) - y) / hx / hy +
           solution[_global(1)] * (x - _x(0)) * (_y(2) - y) / hx / hy +
           solution[_global(2)] * (_x(1) - x) * (y - _y(0)) / hx / hy +
           solution[_global(3)] * (x - _x(0)) * (y - _y(0)) / hx / hy;
}

void FemSolver::set_file(FileType ftype, const std::string & fname)
{
    files[ftype] = fname;
    if (ftype == FT_WORKDIR) files[ftype] = correct_path(fname);
}

std::string FemSolver::correct_path(const std::string & path)
{
    std::string res = path;
    size_t end = res.size() - 1;
    if (res[end] != '/') res.push_back('/');
    return res;
}

inline const std::string FemSolver::file(FileType ftype)
{
    return files[FT_WORKDIR] + files[ftype];
}
