#include "Polymer.h"

using namespace std;

Polymer::Polymer( std::vector< GridSite > gss ) : GridProps( gss[0] ) {
    int maxx = -INT_MAX, maxy = -INT_MAX;
    int minx =  INT_MAX, miny =  INT_MAX;

    atoms = new vector< GridSite >();
    neighbors = new vector< GridSite >();

    for( int i=0; i < gss.size(); i++ ) {
        atoms->push_back( gss[i] );
        maxx = max( maxx, gss[i].X );
        maxy = min( maxy, gss[i].Y );
        minx = max( minx, gss[i].X );
        miny = min( miny, gss[i].Y );
    }
    maxx = max( abs( maxx ), abs( minx ) );
    maxy = max( abs( maxy ), abs( miny ) );

    for( int x= -maxx-1; x <= maxx+1; x++ ) {
        for( int y= -maxy-1; y <= maxy+1; y++ ) {
            GridSite delta( x, y, *this );
            if( !compatiblePosition_lazy( delta ) ) {
                neighbors->push_back( delta );
            }
        }
    }
}

Polymer::~Polymer() {
    atoms->clear();
    delete atoms;
    neighbors->clear();
    delete neighbors;
}

bool Polymer::compatiblePosition_lazy( GridSite delta ) {
    for( int i=0; i < atoms->size(); i++ ) {
        GridSite moved = atoms->at(i) + delta;
        for( int j=0; j < atoms->size(); j++ ) {
            if( moved == atoms->at(i) )
                return false;
        }
    }
    return true;
}

bool Polymer::canStay( Grid<int>& grid, int position ) {
    for( int i=0; i < atoms->size(); i++ ) {
        if( grid[
                ( /*x*/ atoms->at(i).X + position / atoms->at(i).d2 ) % atoms->at(i).d1 * atoms->at(i).d2 + 
                ( /*y*/ atoms->at(i).Y + position ) % atoms->at(i).d2
            ] != GridSite::Free )
            return false;
    }
    return true;
}

void Polymer::depositAndClean( Grid<int>& thegrid, AdvVector& sites, const int position ) const {
    GridSite pos( position, thegrid );
    for( int i=0; i < atoms->size(); i++ ) {
        sites.remove( ( atoms->at(i) + pos ).I() );
        thegrid[ atoms->at(i) + pos ] = GridSite::Atom;
    }
    for( int i=0; i < neighbors->size(); i++ ) {
        sites.remove( ( neighbors->at(i) + pos ).I() );
    }
}

void Polymer::clean( const GridProps& gp, AdvVector& sites, const int position ) const {
    GridSite pos( position, gp );
    for( int i=0; i < atoms->size(); i++ ) {
        sites.remove( ( atoms->at(i) + pos ).I() );
    }
    for( int i=0; i < neighbors->size(); i++ ) {
        sites.remove( ( neighbors->at(i) + pos ).I() );
    }
}

Polymers::Polymers( std::string keyname ) : keyname(keyname), N(0) {}

Polymers::~Polymers() {
    variants.clear();
}

void Polymers::addVariant( Polymer* p ) {
    variants.push_back( p );
    N++;
}

Polymer* Polymers::operator[]( int i ) {
    return variants[i];
}


Polymers* StdPolymers::Dimers::create( GridProps& gp ) {
    Polymers* ps = new Polymers("Dimers");
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite( 0, 1, gp)
    } ) );
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite( 1, 0, gp)
    } ) );
    return ps;
}
Polymers* StdPolymers::LinearTrimers::create( GridProps& gp ) {
    Polymers* ps = new Polymers("Linear trimers");
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite( 0, 1, gp),
        GridSite( 0, 2, gp)
    } ) );
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite( 1, 0, gp),
        GridSite( 2, 0, gp)
    } ) );
    return ps;
}

Polymers* StdPolymers::Trimers::create( GridProps& gp ) {
    Polymers* ps = new Polymers("Trimers");
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite( 0, 1, gp),
        GridSite( 0, 2, gp)
    } ) );
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite( 1, 0, gp),
        GridSite( 2, 0, gp)
    } ) );
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite( 1, 0, gp),
        GridSite( 1, 1, gp)
    } ) );
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite( 0, 1, gp),
        GridSite( 1, 1, gp)
    } ) );
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite(-1, 0, gp),
        GridSite( 0,-1, gp)
    } ) );
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite( 1, 0, gp),
        GridSite( 0, 1, gp)
    } ) );
    return ps;
}

Polymers* StdPolymers::Squared::create( GridProps& gp ) {
    Polymers* ps = new Polymers("Squared");
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite( 0, 1, gp),
        GridSite( 1, 0, gp),
        GridSite( 1, 1, gp)
    } ) );
    return ps;
}

Polymers* StdPolymers::LinearPentamers::create( GridProps& gp ) {
    Polymers* ps = new Polymers("Linear pentamers");
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite( 0, 1, gp),
        GridSite( 0, 2, gp),
        GridSite( 0, 3, gp),
        GridSite( 0, 4, gp)
    } ) );
    ps->addVariant( new Polymer( std::vector< GridSite >{
        GridSite( 0, 0, gp),
        GridSite( 1, 0, gp),
        GridSite( 2, 0, gp),
        GridSite( 3, 0, gp),
        GridSite( 4, 0, gp)
    } ) );
    return ps;
}

std::vector<PolymersFactory*> PolymersFactory::StdPolymers = {
    new StdPolymers::Dimers(),
    new StdPolymers::LinearTrimers(),
    new StdPolymers::Trimers(),
    new StdPolymers::Squared(),
    new StdPolymers::LinearPentamers()
};