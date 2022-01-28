#include "GridFiller.h"
#include "Grid.h"

using namespace std;

void GridFiller::iid(Grid<double>& g) {
    static minstd_rand engine(random_device{}());
    engine.seed( time(NULL) );
    normal_distribution<double> distribution(0.0,1.0);
    for(int i=0; i < g.d1 * g.d2; i++ )
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

int GridFiller::fillWithPolymers(Grid<int>& tofill, Polymers& polys) {

    AdvVector variants( polys.N );
    vector< AdvVector* > sites;
    for( int i=0; i < polys.N; i++ )
        sites.push_back( new AdvVector( tofill.imax() ) );

    // Remove not-free sites
    for( int v=0; v < polys.N; v++ ) {
        for( int i=0; i < tofill.imax(); i++ )
            if( !tofill[i] == GridSite::Free )
                sites[v]->remove(i);

        if( sites[v]->empty() ) variants.remove( v );
    }

    int dep_atoms = 0;
    int var, site;

    // Deposition:
    while( !variants.empty() ) {
        // Select a random variant and ensure sites are available
        var = variants.rnd();

        // Select a random site
        site = sites[var]->rnd();

        // If the site can host the polymer
        if( polys[var]->canStay( tofill, site ) ) {

            // Deposit and remove occupied sites
            for( int i=0; i < polys.N; i++ ) {
                if( i == var )
                    polys[i]->depositAndClean( tofill, *sites[i], site );
                else
                    polys[i]->clean( tofill, *sites[i], site );

                // If this variant cannot be more deposited, remove
                if( sites[i]->empty() )
                    variants.remove( i );
            }

            dep_atoms += polys[var]->atoms->size();
            
        } else {
            
            // Else remove the site
            sites[var]->remove( site );
            if( sites[var]->empty() )
                variants.remove( var );
        }
    }

    // Clean
    for( int i=0; i < polys.N; i++ )
        delete sites[i];

    return dep_atoms;
}