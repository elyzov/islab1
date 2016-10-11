#include "stopwatch.h"
#include <time.h>

Stopwatch::Stopwatch()
{
}

void Stopwatch::start(const char * name)
{
    watchers[name] = clock();
}

real Stopwatch::stop(const char * name)
{
    return (real)(clock() - watchers[name]) / CLOCKS_PER_SEC;
}
