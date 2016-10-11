#ifndef CUBICSPLINEERROR_H
#define CUBICSPLINEERROR_H

#include <iostream>

class CubicSpline_Error
{
public:
    CubicSpline_Error() {}
    virtual void debug_print() const
    {
        std::cerr << "Unknown error" << std::endl;
    }
};

class CubicSpline_OpenFile_Error : public CubicSpline_Error
{
    std::string fname;
public:
    CubicSpline_OpenFile_Error(const std::string & file) : fname(file) {}
    virtual void debug_print() const
    {
        std::cerr << "Can't open file " << fname << std::endl;
    }
};

class CubicSpline_Build_Error : public CubicSpline_Error
{
    std::string msg;
public:
    CubicSpline_Build_Error(const char * message) : msg(message) {}
    virtual void debug_print() const
    {
        std::cerr << "Build exception " << msg << std::endl;
    }
};

#endif // CUBICSPLINEERROR_H
