#include "FourierCoupledGrids.h"

using namespace std;

FourierCoupledGrids::FourierCoupledGrids(int d1) : FourierCoupledGrids(d1,d1) {};

FourierCoupledGrids::FourierCoupledGrids(int d1, int d2) : h(d1,d2), f(d1,d2) {
    direct_plan = fftw_plan_r2r_2d(d1,d2,h.u->data(),f.u->data(), FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
    reverse_plan= fftw_plan_r2r_2d(d1,d2,f.u->data(),h.u->data(), FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
}

void FourierCoupledGrids::fourier_transform(double reverse) {
    if(reverse) {
        fftw_execute( reverse_plan );
        h.normalize();
    } else {
        fftw_execute( direct_plan );
        f.normalize();
    }
}

void FourierCoupledGrids::multiply_fft_old(double gamma) {
    f.u->at(0) = 0;
    pair< int, int > xy;
    double q;
    for(int i=1; i< f.d1 * f.d2; i++) {
        xy = f._xy(i);
        q = sqrt( xy.first * xy.first + xy.second * xy.second );
        f.u->at(i) = f.u->at(i) * pow( q, (gamma-1)/2 );
    }
    f.normalize();
}

void FourierCoupledGrids::multiply_fft_new(double gamma) {
    pair< int, int > xy;
    double q,S;
    double beta = ( gamma - 2 ) / 2;
    double pref = 2 * M_PI / tgamma( beta + 1 );

    for(int i=1; i< f.d1 * f.d2; i++) {
        q = f._abs( _q(i) );
        S = pref * pow( q / 2 , beta ) * cyl_bessel_k( abs(beta), q );
        f.u->at(i) = f.u->at(i) * sqrt( S );
    }
    f.normalize();
}

const std::pair< double, double > FourierCoupledGrids::_q( std::pair< int, int > xy ) const {
    return make_pair(
        2.0 * M_PI / h.d1 * xy.first,
        2.0 * M_PI / h.d2 * xy.second
    );
}

const pair< double, double > FourierCoupledGrids::_q( int i ) const {
    return _q( h._xy(i) );
}