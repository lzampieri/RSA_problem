#include "GridFiller.h"
#include "Grid.h"

using namespace std;

void GridFiller::iid(Grid<double>& g) {
    static minstd_rand engine(random_device{}());
    static bool initialized = false;
    if( !initialized ) {
        initialized = true;
        engine.seed( time(NULL) + hash<thread::id>{}(this_thread::get_id()) );
    }
    static normal_distribution<double> distribution(0.0,1.0);
    for(int i=0; i < g.imax(); i++ )
        g.u->at(i) = distribution(engine);
}

void GridFiller::coscos(Grid<double>& g, double fact1, double fact2) {
    for(int i=0; i < g.d1 * g.d2; i++ ) {
        auto xy = g._xy(i);
        g.u->at(i) = cos( fact1 * xy.X ) * cos( fact2 * xy.Y );
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
        toclean[i] = cleanto;
}

void GridFiller::ranked_insertion(Grid<int>& tofill,const Grid<double>& ranks,const int count) {
    double threshold;
    vector<double> p = *(ranks.u);
    nth_element(p.begin(),p.begin()+count,p.end());
    threshold = p[count];

    // Correct the fact that, if p.begin()+count == p.end(), the function
    // nth_element have no effect
    if( count == p.size() ) {
        threshold = *max_element( p.begin(), p.end() ) + 1;
    }

    int c = 0;    
    for(int i=0; i < ranks.imax(); i++){
        tofill[i] = ranks[i] < threshold ? GridSite::Defect : GridSite::Free;
        c += ( tofill[i] == GridSite::Defect );
    }
    assert( c == count );
}