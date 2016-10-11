#ifndef FEMERROR_H
#define FEMERROR_H

#include <iostream>

class Fem_Error
{
public:
    Fem_Error() {}
    virtual void debug_print() const
    {
        std::cerr << "Unknown error" << std::endl;
    }
};

class Fem_IO_Error : public Fem_Error
{
    std::string msg;
public:
    Fem_IO_Error(const std::string & message) : msg(message) {}
    virtual void debug_print() const
    {
        std::cerr << "Input-Output error: " << msg << std::endl;
    }
};

class Fem_OpenFile_Error : public Fem_Error
{
    std::string fname;
public:
    Fem_OpenFile_Error(const std::string & file) : fname(file) {}
    virtual void debug_print() const
    {
        std::cerr << "Can't open file " << fname << std::endl;
    }
};

class Fem_Grid_Error : public Fem_Error
{
    std::string msg;
public:
    Fem_Grid_Error(const std::string & message) : msg(message) {}
    virtual void debug_print() const
    {
        std::cerr << "Incorrect grid: " << msg << std::endl;
    }
};

class Fem_NVTR_Error : public Fem_Error
{
    std::string msg;
public:
    Fem_NVTR_Error(const std::string & message) : msg(message) {}
    virtual void debug_print() const
    {
        std::cerr << "Error in NVTR: " << msg << std::endl;
    }
};

class Fem_Solver_Error : public Fem_Error
{
    std::string msg;
public:
    Fem_Solver_Error(const std::string & message) : msg(message) {}
    virtual void debug_print() const
    {
        std::cerr << "Error in FemSolver: " << msg << std::endl;
    }
};

#endif // FEMERROR_H
