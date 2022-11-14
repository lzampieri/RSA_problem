#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include "Replicator.h"
#include "Polymer.h"
#include "Grid.h"
#include "NewCF.h"
#include <chrono>
#include "date.h"
#include "AutoScanner.h"

using namespace std;

int main(int argc, char *argv[]) {

    AutoScanner as;
    if( argc > 1 ) {
        string filename = argv[1];
        cout<<"Extracting data from "<<filename<<endl;
        as.loadFromTxt( filename );
    }

    string folder = as.populate();
    auto start_time = chrono::system_clock::now();

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

    auto time_diff = chrono::system_clock::now() - start_time;
    log << "=== Total time: " << date::make_time( time_diff ) << " ===\n";
    cout<< "=== Total time: " << date::make_time( time_diff ) << " ===\n";

}
