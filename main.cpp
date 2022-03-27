#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include "Replicator.h"
#include "Polymer.h"
#include "Grid.h"
#include <chrono>
#include "date.h"
#include "AutoScanner.h"

using namespace std;

int main(int argc, char *argv[]) {

    AutoScanner as;
    if( argc < 2 ) {
        cout<<"Scan params not provided. Proceeding with internal ones."<<endl;
        
        as.chunk_size = 256;
        as.tolerance = 1e-3;
        as.sides = { 256, 512 };
        as.gammas = { 0.5 };
        as.qs = { 0.5 };
        as.ps = {
            new StdPolymers::Trimers()
            };
        as.percolation = true;
        as.draw = false;
        as.verbose = true;
        as.n_threads = 32;
        
    } else {
        string filename = argv[1];
        cout<<"Extracting data from "<<filename<<endl;
        as.loadFromTxt( filename );
    }

    string folder = as.populate();

    ofstream log( folder + "/log.txt", ios_base::app );
    log << "=== " << date::format("%Y%m%d %H:%M", chrono::system_clock::now()) << " ===\n";
    if( as.verbose )
        cout<< "=== " << date::format("%Y%m%d %H:%M", chrono::system_clock::now()) << " ===\n";
    log << ReplicatorParams::header() << endl;
    if( as.verbose )
        cout << ReplicatorParams::header() <<"\n";
    
    int cont = 0;
    for( ReplicatorParams rp : as.rps ) {
        Replicator r( rp );
        log << r.params.to_string() << "(" << ++cont <<" / "<< as.rps.size() << ")\t";
        if( as.verbose )
            cout<< r.params.to_string() << "(" << cont <<" / "<< as.rps.size() << ")\t" << flush;
        r.run();
        log << " Replicas: " << r.total_runned() << endl;
        if( as.verbose )
            cout<< " Replicas: " << r.total_runned() << endl;
        r.save_data();
    }

    system("PAUSE");
}