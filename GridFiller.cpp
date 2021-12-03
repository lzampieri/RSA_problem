#include "GridFiller.h"
#include "Grid.h"

void GridFiller::iid(Grid& g) {
    srand(time(NULL));
    for(int i=0; i < g.d1 * g.d2; i++ )
        g.h->at(i) = ((double)rand()) / RAND_MAX;
}

void GridFiller::coscos(Grid& g, double fact1, double fact2) {
    for(int i=0; i < g.d1 * g.d2; i++ ) {
        auto xy = g._xy(i);
        g.h->at(i) = cos( fact1 * xy.first ) * cos( fact2 * xy.second );
    }
}

void GridFiller::square(Grid& g, int side) {
    if( side < g.d1 && side < g.d2 )
        for( int i= ( g.d1 - side ) / 2; i < ( g.d1 + side) / 2 + 1; i++ ) 
            for( int j= ( g.d2 - side ) / 2; j < ( g.d2 + side) / 2 + 1; j++ )
                g.h->at( i * g.d1 + j ) = 1;
}