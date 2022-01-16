#include <cstdlib>
#include <iostream>
#include <fstream>
#include "Replicator.h"
#include "Polymer.h"
#include "Grid.h"
#include <chrono>
#include "date.h"

using namespace std;

int main() {
    int side = 1024;
    int n_replies = 100;
    bool draw = false;
    GridProps gp( side );
    Polymers* dim = StdPolymers::Dimers( gp );
    Polymers* tri = StdPolymers::LinearTrimers( gp );
    ReplicatorParams rp[] = {
        // Size DefectsFracs Gamma NReplies CorrRange Polymers Draw
        ReplicatorParams( side, 0.2, 0.2, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.2, 0.6, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.2, 1.0, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.2, 1.4, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.3, 0.2, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.3, 0.6, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.3, 1.0, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.3, 1.4, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.4, 0.2, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.4, 0.6, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.4, 1.0, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.4, 1.4, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.5, 0.2, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.5, 0.6, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.5, 1.0, n_replies, -1, dim, draw ),
        ReplicatorParams( side, 0.5, 1.4, n_replies, -1, dim, draw ),

        ReplicatorParams( side, 0.2, 0.2, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.2, 0.6, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.2, 1.0, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.2, 1.4, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.3, 0.2, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.3, 0.6, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.3, 1.0, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.3, 1.4, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.4, 0.2, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.4, 0.6, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.4, 1.0, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.4, 1.4, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.5, 0.2, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.5, 0.6, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.5, 1.0, n_replies, -1, tri, draw ),
        ReplicatorParams( side, 0.5, 1.4, n_replies, -1, tri, draw )
    };

    ofstream log( "log.txt", ios_base::app );
    log << "=== " << date::format("%Y%m%d %H:%M", chrono::system_clock::now()) << " ===\n";
    log << "PATHNAME\t | SIDE | DEFECTS_FRAC | GAMMA | N_REPLIES | CORR_RANGE | POLYMERS\n";
    for( ReplicatorParams r : rp )
        log << Replicator::run_replicator( r ) << endl;

    system("PAUSE");
}