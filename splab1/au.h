#ifndef AUX_H
#define AUX_H

#include "typedef.h"
#include <ostream>

real norm(const rvector & vect);
std::ostream & operator << (std::ostream & out, const rvector & v);

#endif // AUX_H