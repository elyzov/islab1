#include "randomvariable.h"
#include <time.h>
#include <math.h>

RandomVariable::RandomVariable()
{
    eng.seed(time(NULL));
}

real RandomVariable::normalDistribution(real mean, real stddev)
{
    NormalDistribution rv(mean, stddev);
    rv.reset();
    return rv(eng);
}

rvector RandomVariable::normalDistribution(real mean, real stddev, uint count)
{
    rvector res(count);
    NormalDistribution rv(mean, stddev);
    for (auto & r : res) {r = rv(eng); rv.reset();}
    return (rvector&&)res;
}

real RandomVariable::sampleMean(const rvector &vect)
{
    real mean = 0;
    for (auto & val : vect) mean += val;
    return mean / vect.size();
}

real RandomVariable::sampleVariance(const rvector &vect)
{
    real mean = sampleMean(vect);
    real variance = 0;
    for (auto & val : vect) variance += pow(mean - val, 2);
    return variance / (vect.size() - 1);
}

void RandomVariable::setNormalError(const rvector &rate, rvector &vect)
{
    real var = sampleVariance(vect);
    for (uint i = 0; i < vect.size(); ++i)
        vect[i] += normalDistribution(0, rate.at(i) * var);
}

void RandomVariable::setNormalError(real rate, rvector &vect)
{
    real var = sampleVariance(vect);
    rvector nerr = normalDistribution(0, var * rate, vect.size());
    for (uint i = 0; i < vect.size(); ++i)
        vect[i] += nerr[i];
}
