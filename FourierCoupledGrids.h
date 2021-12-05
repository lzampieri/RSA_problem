#ifndef FOURIERCOUPLEDGRIDS_H
#define FOURIERCOUPLEDGRIDS_H

#include <cstdlib>
#include "Grid.h"
#include "fftw3.h"

class FourierCoupledGrids {
protected:
    fftw_plan direct_plan, reverse_plan;

public:
    Grid<double> h;
    Grid<double> f;
    
    FourierCoupledGrids(int d1);
    FourierCoupledGrids(int d1, int d2);

    void fourier_transform(double reverse = false);

    void multiply_fft_old(double gamma);
    void multiply_fft_new(double gamma);

    const std::pair< double, double > _q( std::pair< int, int > xy ) const;
    const std::pair< double, double > _q( int i ) const;

};

#endif