#ifndef POLYMER_H
#define POLYMER_H

#include <vector>
#include <random>
#include "Grid.h"
#include "AdvVector.h"

class Polymer : public GridProps {
public:
    std::vector< int >* atoms;
    std::vector< int >* neighbors;

private:
    bool compatiblePosition_lazy( GridSite delta );
public:
    Polymer( std::vector< GridSite > gss );

    bool canStay( AdvVector& sites, int position );
    void depositAndClean( Grid<int>& thegrid, AdvVector& sites, int position );
};

class Polymers {
public:
    std::vector< Polymer* > variants;
    int N;

    Polymers();
    void addVariant( Polymer* p );
    Polymer* operator[]( int i );
};

namespace StdPolymers {

    Polymers LinearTrimer( GridProps& gp ) {
        Polymers ps;
        ps.addVariant( new Polymer( vector< GridSite >{
            GridSite( 0, 0, gp),
            GridSite( 0, 1, gp),
            GridSite( 0, 2, gp)
        } ) );
        ps.addVariant( new Polymer( vector< GridSite >{
            GridSite( 0, 0, gp),
            GridSite( 1, 0, gp),
            GridSite( 2, 0, gp)
        } ) );
        return ps;
    }

    Polymers Dimers( GridProps& gp ) {
        Polymers ps;
        ps.addVariant( new Polymer( vector< GridSite >{
            GridSite( 0, 0, gp),
            GridSite( 0, 1, gp)
        } ) );
        ps.addVariant( new Polymer( vector< GridSite >{
            GridSite( 0, 0, gp),
            GridSite( 1, 0, gp)
        } ) );
        return ps;
    }

    Polymers Squared( GridProps& gp ) {
        Polymers ps;
        ps.addVariant( new Polymer( vector< GridSite >{
            GridSite( 0, 0, gp),
            GridSite( 0, 1, gp),
            GridSite( 1, 0, gp),
            GridSite( 1, 1, gp)
        } ) );
        return ps;
    }

}

#endif