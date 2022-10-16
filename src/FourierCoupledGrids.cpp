#include "FourierCoupledGrids.h"

using namespace std;

FourierCoupledGrids::FourierCoupledGrids(int d1) : FourierCoupledGrids(d1,d1) {};

FourierCoupledGrids::FourierCoupledGrids(int d1, int d2) : h(d1,d2), f(d1,d2) {
    direct_plan = fftw_plan_r2r_2d(d1,d2,h.u->data(),f.u->data(), FFTW_REDFT00, FFTW_REDFT00, FFTW_MEASURE);
    reverse_plan= fftw_plan_r2r_2d(d1,d2,f.u->data(),h.u->data(), FFTW_REDFT00, FFTW_REDFT00, FFTW_MEASURE);
}

void FourierCoupledGrids::fourier_transform(double reverse) {
    if(reverse) {
        fftw_execute( reverse_plan );
        // h.gaussian_center_and_normalize();
        h.normalize();
    } else {
        fftw_execute( direct_plan );
        f.normalize();
    }
}

void FourierCoupledGrids::multiply_fft_old(double gamma) {
    f.u->at(0) = 0;
    GridSite xy( 0,0,f );
    double q;
    for(int i=1; i< f.d1 * f.d2; i++) {
        q = abs( _q(i) );
        f.u->at(i) = f.u->at(i) * pow( q, (gamma-1)/2 );
    }
    f.normalize();
}

void FourierCoupledGrids::multiply_fft_new(double gamma) {
    pair< int, int > xy;
    double q,S;
    double beta = ( gamma - 2 ) / 2;
    double pref = tgamma( beta + 1 );

    for(int i=0; i< f.d1 * f.d2; i++) {
        q = abs( _q(i) );
        S = pref * pow( q / 2 , beta ) * cyl_bessel_k( abs(beta), q );
        f.u->at(i) = f.u->at(i) * sqrt( S );
    }
    f.normalize();
}

const std::pair< double, double > FourierCoupledGrids::_q( const GridSite& xy ) const {
    double x = xy.X;
    double y = xy.Y;
    if( abs( x ) < 0.05 ) {
        x = 0.05;
    }
    if( abs( y ) < 0.05 ) {
        y = 0.05;
    }
    return make_pair(
        M_PI / ( h.d1 * 2 - 1 ) * x, //( xy.d1 / 2 - abs( xy.X - xy.d1 / 2 ) ),
        M_PI / ( h.d2 * 2 - 1 ) * y  //( xy.d2 / 2 - abs( xy.Y - xy.d2 / 2 ) )
    );
}

const pair< double, double > FourierCoupledGrids::_q( int i ) const {
    return _q( h._xy(i) );
}