#ifndef THEGRID_H
#define THEGRID_H

#define _HAS_STD_BYTE 0

#include <cstdlib>
#include <cmath>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include "fftw3.h"
#include "CorrFuncCalcolator.h"
#include "GridFiller.h"

#define M_PI 3.14159265

class Grid {
protected:
    const int d1, d2;
    std::vector<double>* h;
    std::vector<double>* f;
    fftw_plan direct_plan, reverse_plan;

public:
    Grid(int d1, int d2);
    Grid(int d1);

    void fourier_transform(double reverse = false);

    void multiply_fft_old(double gamma);
    void multiply_fft_new(double gamma);

    void normalize(double fourier = false);

    void print_data(const char* filename, double fourier = false) const;

    // Access
    double operator()(const int i) const; // i \in [0,d1*d2[
    double operator()(const std::pair<int,int> xy, double fourier = false) const; // x \in [0,d1[, y \in [0,d2[
    double operator()(const int x, const int y, double fourier = false) const; // x \in [0,d1[, y \in [0,d2[

protected:
    // PBC
    const int _pbc_x( int x ) const;
    const int _pbc_y( int y ) const;
    const int _pbc_i( int i ) const;
    const std::pair< int, int > _pbc( const std::pair< int, int > xy ) const;
    const std::pair< double, double > _q( std::pair< int, int > xy ) const;
    const std::pair< double, double > _q( int i ) const;
    const double _abs( std::pair< double, double > xy ) const;

    // Coordinates converting
    const double _i(const int x, const int y) const; // x \in [0,d1[, y \in [0,d2[
    const std::pair< int, int > _xy(const int i) const;  // i \in [0,d1*d2[

public:
    // Distances
    const double d(int i1, int i2) const; // i1,i2 \in [0,d1*d2[
    const double d(int x1, int y1, int x2, int y2) const; // x1,x2 \in [0,d1[; y1,y2 \in [0,d2[
    const double d( std::pair< int, int > xy1, std::pair< int, int > xy2 ) const; // x1,x2 \in [0,d1[; y1,y2 \in [0,d2[

    // Friend classes
    friend class CorrFuncCalcolator;
    friend class GridFiller;
};

#endif