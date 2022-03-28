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

GridFiller_Polymers::GridFiller_Polymers( int npolys, int nsites ) {
    variants = new AdvVector( npolys );
    sites = new vector< AdvVector* >();
    for( int i=0; i < npolys; i++ )
        sites->push_back( new AdvVector( nsites ) );
}

GridFiller_Polymers::~GridFiller_Polymers() {
    delete variants;
    for( int i=0; i < sites->size(); i++ )
        delete sites->at(i);
    delete sites;
}

double GridFiller_Polymers::fill(Grid<int>& tofill, Polymers& polys) {
    // Reset
    variants->reset();

    // Remove not available sites
    for( int v=0; v < variants->size; v++ ) {
        sites->at(v)->reset();
        for( int i=0; i < sites->at(v)->size; i++ )
            if( ! ( polys[v]->canStay( tofill, i ) ) )
                sites->at(v)->remove(i);

        if( sites->at(v)->empty() ) variants->remove( v );
    }

    double dep_atoms = 0;
    int var, site;

    // Deposition:
    while( !variants->empty() ) {
        // Select a random variant and ensure sites are available
        var = variants->rnd();

        // Select a random site
        site = sites->at(var)->rnd();

        // If the site can host the polymer
        if( polys[var]->canStay( tofill, site ) ) {

            // Deposit and remove occupied sites
            for( int i=0; i < polys.N; i++ ) {
                if( i == var )
                    polys[i]->depositAndClean( tofill, *sites->at(i), site );
                else
                    polys[i]->clean( tofill, *sites->at(i), site );

                // If this variant cannot be more deposited, remove
                if( sites->at(i)->empty() )
                    variants->remove( i );
            }

            dep_atoms += polys[var]->atoms->size();
            
        } else {
            
            // Else remove the site
            sites->at(var)->remove( site );
            if( sites->at(var)->empty() )
                variants->remove( var );
        }
    }
    return dep_atoms;
}