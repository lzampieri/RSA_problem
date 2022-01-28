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
    int n_replies;
    CorrFunc::Model* CF_model;
    Polymers* to_deposit;
    bool percolation;
    int n_threads;
    bool draw;
    std::string save_path;

    ReplicatorParams( int side, double defects_frac, double gamma,
                      int n_replies, CorrFunc::Model* CF_model = nullptr,
                      Polymers* to_deposit = nullptr, bool percolation = false,
                      int n_threads = 1, bool draw = false,
                      std::string save_path = "" ) :
                      side(side), defects_frac(defects_frac), gamma(gamma),
                      n_replies(n_replies), CF_model(CF_model),
                      to_deposit(to_deposit), percolation(percolation),
                      n_threads(n_threads), draw(draw),
                      save_path(save_path) {};
    std::string to_string();
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
// Conter
    int replicas_to_run;
// Params
    ReplicatorParams params;
// Correlation function
    std::vector< double >* CF_H_avg;
    std::vector< double >* CF_D_avg;
    void update_CF_averages( const std::vector< CorrFunc::Datapoint >* cfh, const std::vector< CorrFunc::Datapoint >* cfd );
// Filling fraction
    double fillfrac_sum;
    double fillfrac_sum2;
    void update_dep_averages( int occupied_sites );
// Percolation
    double defperc_count;
    double atmperc_count;
    void update_perc_averages( bool def_perc, bool atm_perc );
// Works
    std::mutex* mux;
    std::vector< ReplicatorThread* > ongoing;

public:
    Replicator( ReplicatorParams suggested_params );
    ~Replicator();

    void run();
    std::string save_data();

    friend class ReplicatorThread;
};

#endif