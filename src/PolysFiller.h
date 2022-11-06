#ifndef POLYSFILLER_H
#define POLYSFILLER_H

#include "Grid.h"
#include "Polymer.h"
#include "AdvVector.h"
#include <ctime>
#include <cstdlib>
#include <random>
#include <assert.h>


class PolysFiller {
private:
    Grid<int>* grid;
    Polymers* polys;
    AdvVector* available_polysites;

    inline int idx( const int polyid, const int gridid ) const { return polyid * grid->d1d2 + gridid; }
    inline int idx( const int polyid, const GridSite& gridid ) const { return idx( polyid, gridid.I() ); }
    inline int pix( const int polysite ) const { return polysite / grid->d1d2; }
    inline int gix( const int polysite ) const { return polysite % grid->d1d2; }
    inline GridSite gsx( const int polysite ) const { return GridSite( gix( polysite ), *grid ); }
    void clean_around( const GridSite site );

public:
    PolysFiller( Grid<int>* grid, Polymers* polys );
    ~PolysFiller();
    double fill();
};


#endif