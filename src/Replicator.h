#ifndef REPLICATOR_H
#define REPLICATOR_H

#define MAX_REPS 1e8

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <thread>
#include <mutex>
#include <vector>
#include "date.h"
#include "Grid.h"
#include "FourierGrids.h"
#include "NewCF.h"
#include "NewCFH.h"
#include "GridFiller.h"
#include "PolysFiller.h"
#include "Polymer.h"

struct ReplicatorParams {
    int side;
    double defects_frac;
    double gamma;
    int chunk_size;
    double tolerance;
    NewCF::Model* CF_model;
    Polymers* to_deposit;
    int n_threads;
    bool draw;
    bool verbose;
    std::string save_path;

    ReplicatorParams( int side, double defects_frac, double gamma,
                      int chunk_size = 100, double tolerance = 1e-3, NewCF::Model* CF_model = nullptr,
                      Polymers* to_deposit = nullptr,
                      int n_threads = 1, bool draw = false,
                      bool verbose = false, std::string save_path = "" ) :
                      side(side), defects_frac(defects_frac), gamma(gamma),
                      chunk_size(chunk_size), tolerance(tolerance), CF_model(CF_model),
                      to_deposit(to_deposit),
                      n_threads(n_threads), draw(draw),
                      verbose(verbose), save_path(save_path) {};
    std::string to_string();
    static std::string header();
    inline double size() const {  return side*side; }
};

class Replicator;

class ReplicatorThread {
private:
// Infos
    Replicator* thrower;
    int id;
    std::thread* thrd;
// Grid
    Grid<int> g;
    FourierGrids fcg;
// Correlator
    NewCF::Calculator* CFcalc;
    NewCFH::Calculator* CFcalcH;
// Deposition
    PolysFiller* pfl;

// Worker
    static void thread_worker (ReplicatorThread* data);
public:
    ReplicatorThread(Replicator* thrower, int id);
    ~ReplicatorThread();

    void start();
    void join();
};

class Replicator {
private:
// Results
// Correlation function
    std::vector< double >* CF_avg;
    std::vector< double >* CFH_avg;
    void update_NewCF_averages( const std::vector< double >* cf );
    void update_NewCFH_averages( const std::vector< double >* cf );
// Filling fraction
    std::vector< unsigned int >* fills;
    void update_dep_averages( double occupied_sites );
    double fill_avg( unsigned int threshold = UINT_MAX ) const;
    double fill_std( unsigned int threshold = UINT_MAX ) const;
// Works
    std::mutex* mux;
    std::vector< ReplicatorThread* > ongoing;
// Chunks manager
    int runned_replicas;
    int total_replicas_to_run;
    void addChunk();

public:
    Replicator( ReplicatorParams suggested_params );
    ~Replicator();

// Params
    ReplicatorParams params;
    
    void run();
    void save_data();
    int total_runned();

    friend class ReplicatorThread;
};

#endif