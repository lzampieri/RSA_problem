#ifndef CORRFUNCCALCOLATOR_H
#define CORRFUNCCALCOLATOR_H

#include <cstdlib>
#include <vector>
#include <tuple>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <climits>
#include <cstring>
#include "Grid.h"

// #undef _WIN32_WINNT
// #define _WIN32_WINNT 0x0A00  // Windows 10
// #define SRWLOCK_INIT {0}
// #include "mingw_thread/mingw.thread.h"
// #include "mingw_thread/mingw.mutex.h"
#include <thread>
#include <mutex>

#define CORRFUNC_MAX_THREADS 8

#include <iostream>

typedef std::tuple< double, int, double, double > CorrFunc_RawDatapoint; // d, count, sum, sum2
typedef std::tuple< int, double, double > CorrFunc_Datapoint; // d, expval, std

typedef std::pair< int, int > CorrFunc_Work; // v1,v2

template<class T>
class CorrFuncCalcolator {
private:  
    static void _thread_postman(
        const Grid<T>* const grid, const std::vector< CorrFunc_Work* >* const todolist, int& i,
        std::vector< CorrFunc_RawDatapoint >* const results, std::mutex* m );

    static CorrFunc_RawDatapoint _thread_worker(
        const Grid<T>* const grid, const int v1, const int v2 );

public:
    static std::vector< CorrFunc_Datapoint >* compute_corr_function(const Grid<T>* grid, int max_range = INT_MAX); // x, y, std_y

    static void print_corr(const Grid<T>* grid, const char* filename, int max_range = INT_MAX);
};

template class CorrFuncCalcolator<double>;
template class CorrFuncCalcolator<int>;

#endif