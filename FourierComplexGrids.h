#ifndef FOURIERCOMPLEXGRIDS_H
#define FOURIERCOMPLEXGRIDS_H

#include <cstdlib>
#include "Grid.h"
#include "fftw3.h"
#include "utils.h"

class FourierComplexGrids {
protected:
    fftw_plan direct_plan, reverse_plan;

public:
    Grid<double> h;
    fftw_complex* f;
    int half_n1;
    int half_n2;
    int size_of_f;
    
    FourierComplexGrids(int d1);
    FourierComplexGrids(int d1, int d2);

    void fourier_transform(double reverse = false);

    void multiply_fft_old(double gamma);
    void multiply_fft_new(double gamma);

    const double abs_q( int i ) const;
};

#endif