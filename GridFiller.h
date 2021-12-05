#ifndef GRIDFILLER_H
#define GRIDFILLER_H

#include "Grid.h"
#include <ctime>
#include <cstdlib>

class GridFiller {

public:
    static void iid(Grid<double>& g);
    static void coscos(Grid<double>& g, double fact1, double fact2);
    static void square(Grid<double>& g, int side);

    static void clean(Grid<int>& toclean, int cleanto = 0);
    static void montecarlo(Grid<int>& tofill, Grid<double> accept_probs, int count);
};


#endif