#ifndef FOURIERGRIDS_H
#define FOURIERGRIDS_H

#include <cstdlib>
#include "Grid.h"
#include "GridFiller.h"
#include "fftw3.h"
#include "utils.h"

class FourierGrids {
protected:
    fftw_plan direct_plan, reverse_plan;
    fftw_complex* f;

    double* sqrtS;

    int size_of_f;

public:
    Grid<double> h;
    
    FourierGrids(int d1, double gamma);
    FourierGrids(int d1, int d2, double gamma);
    ~FourierGrids();

    void repopulate();

protected:
    void compute_sqrtS( double gamma );

};

#endif