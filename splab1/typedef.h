#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <vector>
#include <list>
#include <string>
#include <utility>
#include <memory>
#include <iostream>
#include <random>
#include <memory>
#include <intrin.h>

#define __PRINT__( val ) (cout << #val << " = " << val << endl)

typedef double real;
typedef unsigned int uint;
typedef unsigned short ushort;
typedef char byte;
typedef unsigned char ubyte;
typedef unsigned long long ullong;

typedef std::vector<real> rvector;
typedef const rvector crvector;
typedef rvector::iterator irvector;
typedef rvector::const_iterator cirvector;

typedef std::vector<int> intVector;
typedef intVector::iterator iintVector;

typedef std::vector<uint> uintVector;

typedef std::list<real> rlist;
typedef std::list<int> ilist;
typedef ilist::iterator _ilist;

typedef std::vector<ilist> ilvector;

typedef std::vector<intVector> intVector2d;
typedef std::vector<rvector> rvector2d;
typedef rvector2d::iterator irvector2d;

typedef std::vector<std::string> stringList;

typedef std::vector<std::pair<int, real> > intRealVector;

typedef std::tr1::ranlux_base_01 RandomEngine;
typedef std::tr1::shared_ptr<RandomEngine> RandomEngine_ptr;
typedef std::tr1::normal_distribution<real> NormalDistribution;
typedef std::tr1::shared_ptr<NormalDistribution> NormalDistribution_ptr;

#endif // TYPEDEF_H
