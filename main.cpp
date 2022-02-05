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
    int sides[] = { 64, 128, 256, 512, 1024, 2048, 4096 };
    double gammas[] = { 0.2, 0.6, 1.0, 1.4, 1.8 };
    double dfs[] = { 0.05, 0.1, 0.2, 0.4 };
    bool percolation = true;
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
        cout<< r.params.to_string() << "(" << ++cont <<" / d"<< rps.size() << ")" << endl;
        r.run();
    }

    system("PAUSE");
}