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
    int side = 100;
    GridProps gp( side );
    Polymers* ps = StdPolymers::Dimers( gp );
    ReplicatorParams rp[] = {
        // Size DefectsFracs Gamma NReplies CorrRange Polymers
        ReplicatorParams( side, 0.2, 0.2, 100, -1, ps ),
        ReplicatorParams( side, 0.2, 0.6, 100, -1, ps ),
        ReplicatorParams( side, 0.2, 1.0, 100, -1, ps ),
        ReplicatorParams( side, 0.2, 1.4, 100, -1, ps ),
        ReplicatorParams( side, 0.3, 0.2, 100, -1, ps ),
        ReplicatorParams( side, 0.3, 0.6, 100, -1, ps ),
        ReplicatorParams( side, 0.3, 1.0, 100, -1, ps ),
        ReplicatorParams( side, 0.3, 1.4, 100, -1, ps ),
        ReplicatorParams( side, 0.4, 0.2, 100, -1, ps ),
        ReplicatorParams( side, 0.4, 0.6, 100, -1, ps ),
        ReplicatorParams( side, 0.4, 1.0, 100, -1, ps ),
        ReplicatorParams( side, 0.4, 1.4, 100, -1, ps ),
        ReplicatorParams( side, 0.5, 0.2, 100, -1, ps ),
        ReplicatorParams( side, 0.5, 0.6, 100, -1, ps ),
        ReplicatorParams( side, 0.5, 1.0, 100, -1, ps ),
        ReplicatorParams( side, 0.5, 1.4, 100, -1, ps )
    };

    ofstream log( "log.txt", ios_base::app );
    log << "=== " << date::format("%Y%m%d %H:%M", chrono::system_clock::now()) << " ===\n";
    log << "PATHNAME\t | SIDE | DEFECTS_FRAC | GAMMA | N_REPLIES | CORR_RANGE | POLYMERS\n";
    for( ReplicatorParams r : rp )
        log << Replicator::run_replicator( r ) << endl;

    system("PAUSE");
}