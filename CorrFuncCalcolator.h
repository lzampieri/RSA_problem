#ifndef CORRFUNCCALCOLATOR_H
#define CORRFUNCCALCOLATOR_H

#include <cstdlib>
#include <vector>
#include <queue>
#include <tuple>
#include <cmath>
// #include <pthread.h>
#include <fstream>

#undef _WIN32_WINNT
#define _WIN32_WINNT 0x0A00  // Windows 10
#include "mingw_thread/mingw.thread.h"
#include "mingw_thread/mingw.mutex.h"

#include <iostream>

using namespace std;

class Grid;

typedef tuple< double, double, double > CorrFunc_datapoint;

class CorrFuncCalcolator {
private:  
    static void _thread_worker(
        const Grid* grid, const int v1, const int v2,
        vector< CorrFunc_datapoint >* results, mutex* m );

public:
    static vector< tuple<double,double, double> >* compute_corr_function(const Grid* grid); // x, y, std_y
    static void print_corr(const Grid* grid, const char* filename);
    static void display_corr(const Grid* grid);

};

#endif