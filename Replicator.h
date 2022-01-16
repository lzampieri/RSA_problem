#ifndef REPLICATOR_H
#define REPLICATOR_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <vector>
#include "date.h"
#include "Grid.h"
#include "FourierCoupledGrids.h"
#include "CorrFuncCalcolator.h"
#include "GridFiller.h"
#include "Polymer.h"

struct ReplicatorParams {
    int side;
    double defects_frac;
    double gamma;
    int n_replies;
    int corr_range;
    Polymers* to_deposit;
    bool draw;
};

class Replicator {
private:
    void run_replica( std::vector< double >* CF_H_avg, std::vector< double >* CF_D_avg, std::vector< int >* pol_dep, std::string draw_file = "" );

protected:
    int side;
    double defects_frac;
    double gamma;
    int n_replies;
    std::string save_path;
    Grid<int> g;
    FourierCoupledGrids fcg;

    // For correlator
    int corr_range;
    CorrFunc::BaseClass<double>* CF_H;
    CorrFunc::BaseClass<int>*    CF_D;

    // For polymers deposition
    Polymers* to_deposit;

    // For drawing
    bool draw;

public:
    Replicator( int side, double defects_frac, double gamma, int n_replies, std::string save_path = "" );
    ~Replicator();

    void enable_correlators( int corr_range = INT_MAX );
    void enable_deposition ( Polymers* to_deposit     );
    void enable_drawing();

    std::string run();
    // Return a string in the form PATHNAME | SIDE | DEFECTS_FRAC | GAMMA | N_REPLIES | CORR_RANGE | POLYMERS

    static std::string run_replicator( ReplicatorParams params );
};

#endif