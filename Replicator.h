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

struct ReplicatorParams {
    int side;
    double defects_frac;
    double gamma;
    int n_replies;
    int corr_range;
};

class Replicator {
private:
    void run_replica( std::vector< double >* CF_H_avg, std::vector< double >* CF_D_avg );

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

public:
    Replicator( int side, double defects_frac, double gamma, int n_replies, std::string save_path = "" );

    void enable_correlators( int corr_range = INT_MAX );

    std::string run();
    // Return a string in the form PATHNAME | SIDE | DEFECTS_FRAC | GAMMA | N_REPLIES |  CORR_RANGE

    static std::string run_replicator( ReplicatorParams params );
};

#endif