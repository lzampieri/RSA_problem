#ifndef UTILS_H
#define UTILS_H

#include <cstdlib>
#include <cmath>

template<class T, class W>
double abs( std::pair<T,W> p ) {
    return sqrt( (double) p.first * p.first + (double) p.second * p.second );
}

#endif