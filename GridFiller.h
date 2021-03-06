#ifndef GRIDFILLER_H
#define GRIDFILLER_H

#include "Grid.h"
#include "Polymer.h"
#include "AdvVector.h"
#include <ctime>
#include <cstdlib>
#include <random>
#include <assert.h>

class GridFiller {
public:
    static void iid(Grid<double>& g);
    static void coscos(Grid<double>& g, double fact1, double fact2);
    static void square(Grid<double>& g, int side);

    static void clean(Grid<int>& toclean, int cleanto = GridSite::Free );
    static void ranked_insertion(Grid<int>& tofill,const Grid<double>& ranks,const int count);

};

class GridFiller_Polymers {
private:
    AdvVector* variants;
    std::vector< AdvVector* >* sites;

public:
    GridFiller_Polymers( int npolys, int nsites );
    ~GridFiller_Polymers();
    double fill(Grid<int>& tofill, Polymers& polys);
};


#endif