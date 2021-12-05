#include "GridFiller.h"
#include "Grid.h"

using namespace std;

void GridFiller::iid(Grid<double>& g) {
    srand(time(NULL));
    for(int i=0; i < g.d1 * g.d2; i++ )
        g.u->at(i) = ((double)rand()) / RAND_MAX;
}

void GridFiller::coscos(Grid<double>& g, double fact1, double fact2) {
    for(int i=0; i < g.d1 * g.d2; i++ ) {
        auto xy = g._xy(i);
        g.u->at(i) = cos( fact1 * xy.first ) * cos( fact2 * xy.second );
    }
}

void GridFiller::square(Grid<double>& g, int side) {
    if( side < g.d1 && side < g.d2 )
        for( int i= ( g.d1 - side ) / 2; i < ( g.d1 + side) / 2 + 1; i++ ) 
            for( int j= ( g.d2 - side ) / 2; j < ( g.d2 + side) / 2 + 1; j++ )
                g.u->at( i * g.d1 + j ) = 1;
}

void GridFiller::clean(Grid<int>& toclean, int cleanto) {
    for( int i=0; i < toclean.d1 * toclean.d2; i++ )
        toclean(i) = cleanto;
}

void GridFiller::montecarlo(Grid<int>& tofill, Grid<double> accept_probs, int count) {
    srand( time( NULL ) );
    while( count ) {
        // Pick a position
        pair<int,int> xy = make_pair(
            rand() % tofill.d1,
            rand() % tofill.d2
        );

        // Check not already filled
        if( tofill(xy) ) continue;

        // Pick a random number
        double accepted = rand();

        // Check if the pick is accepted
        if( accepted > accept_probs( xy ) * RAND_MAX ) {
            // The item is placed!
            tofill(xy) = 1;
            count--;
        }
    }
}