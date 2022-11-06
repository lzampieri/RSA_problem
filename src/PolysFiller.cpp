#include "PolysFiller.h"

using namespace std;

PolysFiller::PolysFiller(Grid<int> *grid, Polymers *polys) : grid(grid), polys(polys)
{
    available_polysites = new AdvVector( grid->d1d2 * polys->N );
}

PolysFiller::~PolysFiller() {
    delete available_polysites;
}

void PolysFiller::clean_around( const GridSite site ) {
    for( int p = 0; p < polys->N; p++ ) {
        for( GridSite g : *polys->at(p)->blocked ) {
            available_polysites->remove( idx( p, site + g ) );
        }
    }
}

double PolysFiller::fill()
{
    // Reset
    available_polysites->reset();

    // Remove not available sites
    for ( int i = 0; i < grid->imax(); i++ )
    {
        if( grid->at(i) != GridSite::Free ) {
            GridSite site( i, *grid );
            clean_around( site );
        }
    }

    double dep_atoms = 0;
    
    int polysite;

    // Deposition:
    while (!available_polysites->empty())
    {
        // Select a random polymer+site
        polysite = available_polysites->rnd();
        Polymer* poly = polys->at( pix( polysite ) );
        GridSite site = gsx( polysite );

        // Assert that the site can host the polymer
        assert( poly->canStay( *grid, site ) );

        // Deposit
        for( const GridSite atom : *poly->atoms ) {
            grid->at( atom + site ) = GridSite::Atom;
            clean_around( atom + site );
        }
        dep_atoms += poly->atoms->size();
    }
    return dep_atoms;
}