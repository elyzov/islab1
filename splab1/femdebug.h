#ifndef FEMDEBUG_H
#define FEMDEBUG_H

#include <iostream>
#include "au.h"
#define LOG(_MSG) FemDebug::inst()(_MSG)

class FemDebug
{
    bool is_debug;

public:
    static FemDebug & inst()
    {
        static FemDebug singleInstance;
        return singleInstance;
    }

    void dubugOn() {is_debug = true;}
    void debugOff() {is_debug = false;}

//    template <class Type>
//    FemDebug & operator << (Type val)
//    {
//        if (is_debug) std::cerr << val;
//        return *this;
//    }

    FemDebug & operator() (std::ostream & out)
    {
        //if (is_debug) std::cerr << out;
        return *this;
    }

private:
    FemDebug();
    FemDebug(const FemDebug &) {}
    const FemDebug& operator = (FemDebug &) {}
};

#endif // FEMDEBUG_H
