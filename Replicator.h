#ifndef REPLICATOR_H
#define REPLICATOR_H

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
#include "FourierCoupledGrids.h"
#include "CorrFuncCalcolator.h"
#include "GridFiller.h"
#include "Polymer.h"
#include "Percolator.h"

struct ReplicatorParams {
    int side;
    double defects_frac;
    double gamma;
    int chunk_size;
    double tolerance;
    CorrFunc::Model* CF_model;
    Polymers* to_deposit;
    bool percolation;
    int n_threads;
    bool draw;
    bool verbose;
    std::string save_path;

    ReplicatorParams( int side, double defects_frac, double gamma,
                      int chunk_size = 100, double tolerance = 1e-3, CorrFunc::Model* CF_model = nullptr,
                      Polymers* to_deposit = nullptr, bool percolation = false,
                      int n_threads = 1, bool draw = false,
                      bool verbose = false, std::string save_path = "" ) :
                      side(side), defects_frac(defects_frac), gamma(gamma),
                      chunk_size(chunk_size), tolerance(tolerance), CF_model(CF_model),
                      to_deposit(to_deposit), percolation(percolation),
                      n_threads(n_threads), draw(draw),
                      verbose(verbose), save_path(save_path) {};
    std::string to_string();
    static std::string header();
    double size() const;
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
    FourierCoupledGrids fcg;
    GridFiller_Polymers gfp;
// Correlator
    CorrFunc::Calculator<double>* CF_H;
    CorrFunc::Calculator<int>*    CF_D;
// Percolation
    Percolator* perc;

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
    std::vector< double >* CF_H_avg;
    std::vector< double >* CF_D_avg;
    void update_CF_averages( const std::vector< CorrFunc::Datapoint >* cfh, const std::vector< CorrFunc::Datapoint >* cfd );
// Filling fraction
    std::vector< unsigned int >* fills;
    void update_dep_averages( double occupied_sites );
    double fill_avg( unsigned int threshold = UINT_MAX ) const;
    double fill_std( unsigned int threshold = UINT_MAX ) const;
    double fill_std_fromfit( unsigned int threshold = UINT_MAX ) const;
// Percolation
    double defperc_count;
    double atmperc_count;
    void update_perc_averages( bool def_perc, bool atm_perc );
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