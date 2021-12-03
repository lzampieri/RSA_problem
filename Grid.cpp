#include "Grid.h"

using namespace std;

Grid::Grid(int d1) : Grid(d1, d1) {};

Grid::Grid(int d1, int d2) : d1(d1), d2(d2) {
    h = new vector<double>(d1*d2, 0);
    f = new vector<double>(d1*d2, 0);
    direct_plan = fftw_plan_r2r_2d(d1,d2,h->data(),f->data(), FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
    reverse_plan= fftw_plan_r2r_2d(d1,d2,f->data(),h->data(), FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
}

void Grid::fourier_transform(double reverse) {
    if(reverse) fftw_execute( reverse_plan );
    else fftw_execute( direct_plan );
    normalize( !reverse );
}

void Grid::multiply_fft_old(double gamma) {
    f->at(0) = 0;
    pair< int, int > xy;
    double q;
    for(int i=1; i<d1*d2; i++) {
        xy = _xy(i);
        q = sqrt( xy.first * xy.first + xy.second * xy.second );
        f->at(i) = f->at(i) * pow( q, (gamma-1)/2 );
    }
    normalize( true );
}

void Grid::multiply_fft_new(double gamma) {
    pair< int, int > xy;
    double q,S;
    double beta = ( gamma - 2 ) / 2;
    double pref = 2 * M_PI / tgamma( beta + 1 );

    for(int i=1; i<d1*d2; i++) {
        q = _abs( _q(i) );
        S = pref * pow( q / 2 , beta ) * cyl_bessel_k( abs(beta), q );
        f->at(i) = f->at(i) * sqrt( S );
    }
    normalize( true );
}

void Grid::normalize(double fourier) {
    vector<double>* v = fourier ? f : h;
    double maximum = v->at(0);
    for(int i=1; i<d1*d2; i++) {
        maximum = max( maximum, v->at(i) );
    }
    for(int i=0; i<d1*d2; i++) {
        v->at(i) /= maximum;
    }
}

void Grid::print_data(const char* filename, double fourier) const {
    ofstream out(filename);
    for(int i=0; i<d1; i++)
        for(int j=0; j<d2; j++)
            out<<i<<'\t'<<j<<'\t'<<operator()(i,j,fourier)<<'\n';
    out.close();
}

// Access

double Grid::operator()(const int i) const {
    return h->at( _pbc_i(i) );
}

double Grid::operator()(const pair<int,int> xy, double fourier) const {
    return operator()( xy.first, xy.second, fourier );
}

double Grid::operator()(const int x, const int y, double fourier) const {
    if( fourier ) return f->at( _i(x,y) );
    return h->at( _i(x,y) );
}

// PBC

const int Grid::_pbc_x( int x ) const {
    while( x < 0 ) x += d1;
    return x % d1;
}

const int Grid::_pbc_y( int y ) const {
    while( y < 0 ) y += d1;
    return y % d2;
}

const int Grid::_pbc_i( int i ) const {
    while( i < 0 ) i += d1*d2;
    return i % (d1*d2);
}

const pair< int, int > Grid::_pbc( const pair< int, int > xy ) const {
    return make_pair( _pbc_x( xy.first ), _pbc_y( xy.second ) );
}

const pair< double, double > Grid::_q( pair< int, int > xy ) const {
    return make_pair(
        2.0 * M_PI / d1 * xy.first,
        2.0 * M_PI / d2 * xy.second
    );
}

const pair< double, double > Grid::_q( int i ) const {
    return _q( _xy(i) );
}

const double Grid::_abs( pair< double, double > xy ) const {
    return sqrt( xy.first * xy.first + xy.second * xy.second );
}

// Coordinates converting

const double Grid::_i(const int x, const int y) const {
    return _pbc_x(x) * d2 + _pbc_y(y);
}

const pair< int, int > Grid::_xy(const int i) const {
    return make_pair< int, int >( _pbc_i(i) / d2, _pbc_i(i) % d2 );   
}

// Distances

const double Grid::d(int i1, int i2) const{
    return d( _xy(i1), _xy(i2) );
}

const double Grid::d(int x1, int y1, int x2, int y2) const{
    int deltaX = _pbc_x( x1 - x2 );
    deltaX = min( deltaX, d1 - deltaX );
    int deltaY = _pbc_x( y1 - y2 );
    deltaY = min( deltaY, d2 - deltaY );
    return sqrt( pow( deltaX, 2 ) + pow( deltaY, 2 ) );
}

const double Grid::d( pair< int, int > xy1, pair< int, int > xy2 ) const {
    return d( xy1.first, xy1.second, xy2.first, xy2.second );
}