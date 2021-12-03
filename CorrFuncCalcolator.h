#ifndef CORRFUNCCALCOLATOR_H
#define CORRFUNCCALCOLATOR_H

#include <cstdlib>
#include <vector>
#include <queue>
#include <tuple>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <climits>

#undef _WIN32_WINNT
#define _WIN32_WINNT 0x0A00  // Windows 10
#include "mingw_thread/mingw.thread.h"
#include "mingw_thread/mingw.mutex.h"

#define CORRFUNC_MAX_THREADS 8

#include <iostream>

class Grid;

typedef std::tuple< double, int, double, double > CorrFunc_RawDatapoint; // d, count, sum, sum2
typedef std::tuple< int, double, double > CorrFunc_Datapoint; // d, expval, std

typedef std::pair< int, int > CorrFunc_Work; // v1,v2

class CorrFuncCalcolator {
private:  
    static void _thread_postman(
        const Grid* grid, std::queue< CorrFunc_Work >* todolist,
        std::vector< CorrFunc_RawDatapoint >* results, std::mutex* m );
    static CorrFunc_RawDatapoint _thread_worker(
        const Grid* grid, const int v1, const int v2 );

public:
    static std::vector< CorrFunc_Datapoint >* compute_corr_function(const Grid* grid, int max_range = INT_MAX); // x, y, std_y
    static void print_corr(const Grid* grid, const char* filename, int max_range = INT_MAX);
    static void display_corr(const Grid* grid, int max_range = INT_MAX);

};

#endif