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
    int n_replies = 1024;
    int sides[] = { 2048 };
    double gammas[] = { 0.4, 0.8, 1.2, 1.6 };
    double dfs[] = { 0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. };
    bool percolation = true;
    bool draw = false;
    int n_threads = 32;

    vector< ReplicatorParams > rps;
    for( int s : sides ) {
        GridProps gp( s );
        // Polymers* pols[] = { StdPolymers::Dimers( gp ), StdPolymers::Trimers( gp ), StdPolymers::LinearTrimers( gp ), StdPolymers::Squared( gp ) };
        Polymers* pols[] = { StdPolymers::Dimers( gp ) };
        for( double g : gammas )
            for( double df : dfs )
                for( Polymers* p : pols )
                    // Size DefectsFracs Gamma NReplies CorrRange Polymers Draw
                    rps.push_back( ReplicatorParams( s, df, g, n_replies, nullptr, p, percolation, n_threads, draw ) );
    }

    ofstream log( "log.txt", ios_base::app );
    log << "=== " << date::format("%Y%m%d %H:%M", chrono::system_clock::now()) << " ===\n";
    cout<< "=== " << date::format("%Y%m%d %H:%M", chrono::system_clock::now()) << " ===\n";
    log << ReplicatorParams::header() <<"\n";
    cout << ReplicatorParams::header() <<"\n";
    
    int cont = 0;
    for( ReplicatorParams rp : rps ) {
        Replicator r( rp );
        log << r.params.to_string() << '\n';
        cout<< r.params.to_string() << "(" << ++cont <<" / "<< rps.size() << ")" << endl;
        r.run();
        r.save_data();
    }

    system("PAUSE");
}