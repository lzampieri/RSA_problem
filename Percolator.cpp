#include "Percolator.h"

using namespace std;

Percolator::Percolator(Grid<int>& grid) : grid(grid) {
    reset();
}

bool Percolator::next(GridSite& i) {
    // Remember than in c++:
    // operator++( i ) -> i += 1; return i-1;
    if( to_check_it == 0 ) return false;
    i = to_check[ --to_check_it ];
    return true;
}

bool Percolator::add(GridSite i) {
    if( inserted[ i.I() ] ) return false;
    inserted[ i.I() ] = true;
    to_check[ to_check_it++ ] = i;
    return true;
}

void Percolator::reset() {
    inserted.resize( grid.imax(), false );
    fill( inserted.begin(), inserted.end(), false );
    to_check.resize( grid.imax(), GridSite(0, grid) );

    // Add the first row of the grid;
    for( int i=0; i < grid.d1; i++ ) {
        to_check[i] = GridSite( i, grid );
        inserted[i] = true;
    }
    to_check_it = grid.d1;
}

bool Percolator::is_percolating(short what) {
    reset();
    GridSite i(0, grid);
    GridSite l( 1, 0, grid ), r( -1, 0, grid ), d( 0, 1, grid );
    while( next(i) ) {
        if( grid[ i ] == what ) {
            if( i.last_row() ) return true;
            add( i + l );
            add( i + r );
            add( i + d );
        }
    }
    return false;
}