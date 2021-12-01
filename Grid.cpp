#include "Grid.h"

// Todo sistemare

Grid::Grid(int d1) : Grid(d1, d1) {};

Grid::Grid(int d1, int d2) : d1(d1), d2(d2) {
    h = new vector<double>(d1*d2, 0);
    f = new vector<double>(d1*d2, 0);
    direct_plan = fftw_plan_r2r_2d(d1,d2,h->data(),f->data(), FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
    reverse_plan= fftw_plan_r2r_2d(d1,d2,f->data(),h->data(), FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
}

void Grid::populate_iid() {
    srand(time(NULL));
    for(int i=0; i < d1*d2; i++ )
        h->at(i) = ((double)rand()) / RAND_MAX;
}

void Grid::populate_coscos(double fact1, double fact2) {
    for(int i=0; i < d1*d2; i++ )
        h->at(i) = cos( fact1 * ( i % d1 ) ) * cos( fact2 * ( i / d1 ) );
}

void Grid::populate_square(int side) {
    for( int i= ( d1 - side ) / 2; i < ( d1 + side) / 2 + 1; i++ ) 
        for( int j= ( d2 - side ) / 2; j < ( d2 + side) / 2 + 1; j++ )
            h->at( i * d1 + j ) = 0.5;
}

void Grid::fourier_transform(double reverse) {
    if(reverse) fftw_execute( reverse_plan );
    else fftw_execute( direct_plan );
}

void Grid::multiply_fft(double gamma) {
    f->at(0) = 0;
    for(int i=1; i<d1*d2; i++)
        f->at(i) = f->at(i) * pow( d(i,0), (gamma-1)/2 );
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

const int Grid::_pbc_x( const int x ) const {
    return x % d1;
}

const int Grid::_pbc_y( const int y ) const {
    return y % d2;
}

const int Grid::_pbc_i( const int i ) const {
    return i % (d1*d2);
}

const pair< int, int > Grid::_pbc( const pair< int, int > xy ) const {
    return make_pair( _pbc_x( xy.first ), _pbc_y( xy.second ) );
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