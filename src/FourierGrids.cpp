#include "FourierGrids.h"

using namespace std;

FourierGrids::FourierGrids(int d1, double gamma) : FourierGrids(d1,d1,gamma) {};

FourierGrids::FourierGrids(int d1, int d2, double gamma) : h(d1,d2) {
    
    // Init g and plans
    size_of_f = d1 * ( d2 / 2 + 1 );
    f = fftw_alloc_complex( size_of_f );
    direct_plan = fftw_plan_dft_r2c_2d( d1, d2, h.u->data(), f, FFTW_MEASURE);
    reverse_plan= fftw_plan_dft_c2r_2d( d1, d2, f, h.u->data(), FFTW_MEASURE);

    // Init sqrtS
    sqrtS = new double[ size_of_f ];
    compute_sqrtS( gamma );
}

FourierGrids::~FourierGrids() {
    fftw_destroy_plan( direct_plan );
    fftw_destroy_plan( reverse_plan );
    fftw_free( f );
    delete[] sqrtS;
}

void FourierGrids::compute_sqrtS( double gamma ) {
    // Create a new grid
    Grid<double> grid( h.d1, h.d2 );
    double dist_x, dist_y, dist;
    for( int x = 0; x < grid.d1; x++ ) {
        for( int y = 0; y < grid.d2; y++ ) {
            dist_x = min( x, grid.d1 - x );
            dist_y = min( y, grid.d2 - y );
            dist = sqrt( dist_x * dist_x + dist_y * dist_y );
            grid( x, y ) = pow( 1.0 + dist * dist, - gamma / 2 );
        }
    }

    // Fourier transform
    fftw_complex* S = fftw_alloc_complex( size_of_f );
    fftw_plan S_plan = fftw_plan_dft_r2c_2d( grid.d1, grid.d2, grid.u->data(), S, FFTW_ESTIMATE );
    fftw_execute( S_plan );
    
    // Compute square root
    for( int i = 0; i < size_of_f; i++ ) {
        sqrtS[i] = sqrt( S[i][0] );
    }

    // clean
    fftw_destroy_plan( S_plan );
    fftw_free( S );
}

void FourierGrids::repopulate() {

    // Fill h with random variables
    GridFiller::iid( h );
    
    // Compute the fourier transform
    fftw_execute( direct_plan );

    // Multiply by fftw
    for( int i = 0; i < size_of_f; i++ ) {
        f[ i ][ 0 ] *= sqrtS[ i ];
        f[ i ][ 1 ] *= sqrtS[ i ];
    }

    // Reverse the fourier transform
    fftw_execute( reverse_plan );
}