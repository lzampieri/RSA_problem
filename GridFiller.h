#ifndef GRIDFILLER_H
#define GRIDFILLER_H

#include "Grid.h"
#include <ctime>
#include <cstdlib>
#include <random>

class GridFiller {
public:
    static void iid(Grid<double>& g);
    static void coscos(Grid<double>& g, double fact1, double fact2);
    static void square(Grid<double>& g, int side);

    static void clean(Grid<int>& toclean, int cleanto = GridSite::Free );
    static void ranked_insertion(Grid<int>& tofill,const Grid<double> ranks,const int count);
};


#endif