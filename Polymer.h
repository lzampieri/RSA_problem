#ifndef POLYMER_H
#define POLYMER_H

#include <vector>
#include <random>
#include <string>
#include "Grid.h"
#include "AdvVector.h"

class Polymer : public GridProps {
public:
    std::vector< GridSite >* atoms;
    std::vector< GridSite >* neighbors;

private:
    bool compatiblePosition_lazy( GridSite delta );
public:
    Polymer( std::vector< GridSite > gss );

    bool canStay( Grid<int>& grid, int position );
    void depositAndClean( Grid<int>& thegrid, AdvVector& sites, int position );
    void clean( const GridProps& gp, AdvVector& sites, int position );
};

class Polymers {
public:
    std::vector< Polymer* > variants;
    int N;
    std::string keyname;

    Polymers( std::string keyname );
    void addVariant( Polymer* p );
    Polymer* operator[]( int i );

    Polymers(const Polymers&) = delete; // Copy costruction forbidden
};

namespace StdPolymers {
    Polymers* LinearTrimers( GridProps& gp );
    Polymers* LinearPentamers( GridProps& gp );
    Polymers* Dimers( GridProps& gp );
    Polymers* Squared( GridProps& gp );
}

#endif