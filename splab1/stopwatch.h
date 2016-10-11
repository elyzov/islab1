#ifndef STOPWATCH_H
#define STOPWATCH_H

#include <map>
#include "typedef.h"


typedef std::map<const char *, ullong> StapwatchTimeMap;

class Stopwatch
{
    StapwatchTimeMap watchers;
public:
    static Stopwatch & inst()
    {
        static Stopwatch singleInstance;
        return singleInstance;
    }

    void start(const char *name);
    real stop(const char *name);

    ullong rdtsc()
    {
		return __rdtsc();
    }

private:
    Stopwatch();
    Stopwatch(const Stopwatch &) {}
    const Stopwatch& operator = (Stopwatch &) {}
};

#endif // STOPWATCH_H
