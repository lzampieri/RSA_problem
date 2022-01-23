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
    int n_replies = 100;
    int sides[] = { 64, 128, 256, 512, 1024, 2048 };
    double gammas[] = { 0.2, 0.6, 1.0, 1.4 };
    double dfs[] = { 0.2, 0.3, 0.4, 0.5 };
    bool draw = false;
    int n_threads = 64;

    vector< ReplicatorParams > rps;
    for( int s : sides ) {
        GridProps gp( s );
        Polymers* pols[] = { StdPolymers::Dimers( gp ), StdPolymers::Trimers( gp ), StdPolymers::LinearTrimers( gp ), StdPolymers::Squared( gp ) };
        for( double g : gammas )
            for( double df : dfs )
                for( Polymers* p : pols )
                    // Size DefectsFracs Gamma NReplies CorrRange Polymers Draw
                    rps.push_back( ReplicatorParams( s, df, g, n_replies, nullptr, p, n_threads, draw ) );
    }

    ofstream log( "log.txt", ios_base::app );
    log << "=== " << date::format("%Y%m%d %H:%M", chrono::system_clock::now()) << " ===\n";
    cout<< "=== " << date::format("%Y%m%d %H:%M", chrono::system_clock::now()) << " ===\n";
    log << "PATHNAME\t | SIDE | DEFECTS_FRAC | GAMMA | N_REPLIES | CORR_RANGE | POLYMERS\n";
    cout<< "[log]        PATHNAME\t | SIDE | DEFECTS_FRAC | GAMMA | N_REPLIES | CORR_RANGE | POLYMERS\n";
    for( ReplicatorParams rp : rps ) {
        Replicator r( rp );
        r.run();
        log << r.save_data() << endl;
    }

    system("PAUSE");
}