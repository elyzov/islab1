#include <iostream>
#include <stdlib.h>
#include "femdatgen.h"
#include "femsolver.h"
#include "typedef.h"
#include "feminverse.h"
#include "denseslae.h"
#include "randomvariable.h"
#include "femresonspline.h"
#include "au.h"

//namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[]) try
{
    string inputFile = "";
    string pointsFile = "points.txt";
    string pathOut = ".";
    string splineFile = "";
    string linear, nonlinear;
    bool withOutput = false;

    for (short i = 1; i < argc; i++)
    {
        string token(argv[i]);
        if (token == "-i")
        {
            if (++i < argc) inputFile = argv[i];
            else throw Fem_IO_Error("не заданно значение для флага -i");
        }
        else if (token == "-r")
        {
            if (++i < argc) splineFile = argv[i];
            else throw Fem_IO_Error("не заданно значение для флага -r");
        }
        else if (token == "-p")
        {
            if (++i < argc) pointsFile = argv[i];
            else throw Fem_IO_Error("не заданно значение для флага -p");
        }
        else if (token == "-o")
        {
            withOutput = true;
            if (++i < argc) pathOut = argv[i];
            else throw Fem_IO_Error("не заданно значение для флага -o");
        }
        else if (token == "-l")
        {
            if (++i < argc) linear = argv[i];
            else throw Fem_IO_Error("не заданно значение для флага -l");
        }
        else if (token == "-n")
        {
            if (++i < argc) nonlinear = argv[i];
            else throw Fem_IO_Error("не заданно значение для флага -n");
        }
        else throw Fem_IO_Error("неизвестный аргумент");
    }

    if (!inputFile.empty())
    {
        if (linear.empty() && nonlinear.empty())
            throw Fem_IO_Error("Не определена обратная задача.");
        FemDatGen gen(inputFile);
        FemSolver dirTask(gen, pointsFile, pathOut);
        if (withOutput)
        {
            gen.dataOutput(pathOut);
            dirTask.solve();
            dirTask.output();
        }
        FemInverse invTask(&dirTask);
        invTask.solveLinear(linear);
        invTask.solveNonLinear(nonlinear);
    }
    else if (!splineFile.empty())
    {
        FemResOnSpline task(splineFile, pointsFile);
        rvector sol;
        task.diffSolInPointsPair(sol);
        cout << "res = " << sol << endl;
    }
    else throw Fem_IO_Error("Не определена задача.");
    exit(0);
}
/*! --------------------- Обработка ошибок --------------------- */
catch (Fem_Error & error)
{
    std::cerr << "FEM Error: ";
    error.debug_print();
}
catch (CubicSpline_Error & error)
{
    std::cerr << "CubicSpline Error: ";
    error.debug_print();
}
catch (DenseSLAE_Error & error)
{
    std::cerr << "DenseSLAE Error: ";
    error.print_debug();
}

catch (std::ifstream::failure &)
{
    std::cerr << "Input error: Incorrect format of file" << std::endl;
}

catch (std::exception & error)
{
    std::cerr << "STL Exception: " << error.what() << std::endl;
}
catch (...)
{
    std::cerr << "Unknown error" << std::endl;
}
