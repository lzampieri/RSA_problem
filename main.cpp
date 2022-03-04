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
    int chunk_size = 1024;
    double tolerance = 1e-5;
    // int sides[] = { 2048 };
    // double gammas[] = { 0.4, 0.8, 1.2, 1.6 };
    // double dfs[] = { 0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. };
    int sides[] = { 64, 128, 256 };
    double gammas[] = { 0.8 };
    double dfs[] = { 0.3 };
    bool percolation = true;
    bool draw = false;
    bool verbose = true;
    int n_threads = 32;

    vector< ReplicatorParams > rps;
    for( int s : sides ) {
        GridProps gp( s );
        // Polymers* pols[] = { StdPolymers::Dimers( gp ), StdPolymers::Trimers( gp ), StdPolymers::LinearTrimers( gp ), StdPolymers::Squared( gp ) };
        Polymers* pols[] = { StdPolymers::Dimers( gp ) };
        for( double g : gammas )
            for( double df : dfs )
                for( Polymers* p : pols )
                    //                               Size DefectsFracs Gamma ChunkSize   Tolerance  CFModel  Polymers Percolation  NThreads,  Draw  Verbose SavePath
                    rps.push_back( ReplicatorParams( s,   df,          g,    chunk_size, tolerance, nullptr, p,       percolation, n_threads, draw, verbose          ) );
    }

    ofstream log( "log.txt", ios_base::app );
    log << "=== " << date::format("%Y%m%d %H:%M", chrono::system_clock::now()) << " ===\n";
    if( verbose )
        cout<< "=== " << date::format("%Y%m%d %H:%M", chrono::system_clock::now()) << " ===\n";
    log << ReplicatorParams::header() <<"\n";
    if( verbose )
        cout << ReplicatorParams::header() <<"\n";
    
    int cont = 0;
    for( ReplicatorParams rp : rps ) {
        Replicator r( rp );
        log << r.params.to_string() << '\n';
        if( verbose )
            cout<< r.params.to_string() << "(" << ++cont <<" / "<< rps.size() << ")" << endl;
        r.run();
        log << "// Chunks: " << r.runned_chunks() << '\n';
        if( verbose )
            cout<< "// Chunks: " << r.runned_chunks() << endl;
        r.save_data();
    }

    system("PAUSE");
}