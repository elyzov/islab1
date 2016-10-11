#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include "typedef.h"

class RandomVariable
{
    RandomEngine eng;

public:
    static RandomVariable & inst()
    {
        static RandomVariable singleInstance;
        return singleInstance;
    }

    real normalDistribution(real mean, real stddev);
    rvector normalDistribution(real mean, real stddev, uint count);
    real sampleMean(const rvector & vect);
    real sampleVariance(const rvector & vect);
    void setNormalError(const rvector & rate, rvector & vect);
    void setNormalError(real rate, rvector & vect);

private:
    RandomVariable();
    RandomVariable(const RandomVariable &) {}
    const RandomVariable& operator = (RandomVariable &) {}
};

#endif // RANDOMVARIABLE_H
