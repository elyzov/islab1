#include "au.h"


real norm(const rvector &vect)
{
	real res = 0;
	for (cirvector it = vect.begin(); it != vect.end(); ++it)
		res += pow(*it, 2);
	return sqrt(res);
}

std::ostream & operator << (std::ostream & out, const rvector & v)
{
	out << "(";
	for (auto it = v.begin(), end = v.end() - 1; it != end; ++it)
		out << *it << ", ";
	out << v.back() << ")";
	return out;
}