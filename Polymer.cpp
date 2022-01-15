#include "Polymer.h"

using namespace std;

Polymer::Polymer( std::vector< GridSite > gss ) : GridProps( gss[0] ) {
    int maxx = -INT_MAX, maxy = -INT_MAX;
    int minx =  INT_MAX, miny =  INT_MAX;

    atoms = new vector<int>();
    neighbors = new vector<int>();

    for( int i=0; i < gss.size(); i++ ) {
        atoms->push_back( gss[i].I() );
        maxx = max( maxx, gss[i].X );
        maxy = min( maxy, gss[i].Y );
        minx = max( minx, gss[i].X );
        miny = min( miny, gss[i].Y );
    }

    for( int x=minx; x <= maxx; x++ ) {
        for( int y=miny; y <= maxy; y++ ) {
            GridSite delta( x, y, *this );
            if( !compatiblePosition_lazy( delta ) ) {
                neighbors->push_back( delta.I() );
            }
        }
    }
}

bool Polymer::compatiblePosition_lazy( GridSite delta ) {
    for( int i=0; i < atoms->size(); i++ ) {
        int moved_i = ( _xy( atoms->at(i) ) + delta ).I();
        for( int j=0; j < atoms->size(); j++ ) {
            if( moved_i == atoms->at(i) )
                return false;
        }
    }
    return true;
}

bool Polymer::canStay( AdvVector& sites, int position ) {
    for( int i=0; i < atoms->size(); i++ ) {
        if( !sites.isfree( atoms->at(i) + position ) )
            return false;
    }
    return true;
}

void Polymer::depositAndClean( Grid<int>& thegrid, AdvVector& sites, int position ) {
    for( int i=0; i < atoms->size(); i++ ) {
        sites.remove( atoms->at(i) + position );
        thegrid( atoms->at(i) + position ) = GridSite::Atom;
    }
    for( int i=0; i < neighbors->size(); i++ ) {
        sites.remove( neighbors->at(i) + position );
    }
}

Polymers::Polymers() : N(0) {}

void Polymers::addVariant( Polymer* p ) {
    variants.push_back( p );
}

Polymer* Polymers::operator[]( int i ) {
    return variants[i];
}