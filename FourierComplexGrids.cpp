#include "FourierComplexGrids.h"

using namespace std;

FourierComplexGrids::FourierComplexGrids(int d1) : FourierComplexGrids(d1,d1) {};

FourierComplexGrids::FourierComplexGrids(int d1, int d2) : h(d1,d2) {
    half_n1 = ( d1 / 2 + 1 );
    half_n2 = ( d2 / 2 + 1 );
    size_of_f = d1 * half_n2;
    f = fftw_alloc_complex( size_of_f );
    direct_plan = fftw_plan_dft_r2c_2d( d1, d2, h.u->data(), f, FFTW_MEASURE);
    reverse_plan= fftw_plan_dft_c2r_2d( d1, d2, f, h.u->data(), FFTW_MEASURE);
}

void FourierComplexGrids::fourier_transform(double reverse) {
    if(reverse) {
        fftw_execute( reverse_plan );
    } else {
        fftw_execute( direct_plan );
    }
}

void FourierComplexGrids::multiply_fft_old(double gamma) {
    f[0][0] = 0;
    f[0][1] = 0;
    double q;
    for(int i=1; i< size_of_f; i++) {
        q = abs_q(i);
        f[i][0] *= pow( q, (gamma-1)/2 );
        f[i][1] *= pow( q, (gamma-1)/2 );
    }
}

void FourierComplexGrids::multiply_fft_new(double gamma) {
    pair< int, int > xy;
    double q,S;
    double beta = ( gamma - 2 ) / 2;
    double pref = tgamma( beta + 1 );

    for(int i=0; i<size_of_f; i++) {
        q = abs_q(i);
        S = pref * pow( q / 2 , beta ) * cyl_bessel_k( abs(beta), q );
        f[i][0] *= sqrt( S );
        f[i][1] *= sqrt( S );
    }
}

const double FourierComplexGrids::abs_q( int i ) const {
    double x = floor( i / half_n2 );
    if( x > half_n1 )
        x = 2 * half_n1 - x;
    double y = i % half_n2;
    if( ( abs( x ) < 0.05 ) && ( abs( y ) < 0.05 ) ) {
        x = 0.1;
        y = 0.1;
    }
    double q_x = 2 * M_PI / half_n1 * x;
    double q_y = 2 * M_PI / half_n2 * y;
    return sqrt( q_x * q_x + q_y * q_y );
}