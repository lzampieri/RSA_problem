#include "AdvVector.h"

using namespace std;

minstd_rand AdvVector::engine = minstd_rand{ std::random_device{}() };

AdvVector::AdvVector( int size ) : size( size ) {
    data = new int[size];
    whereis = new int[size];
    for( int i=0; i < size; i++ ) {
        data[i] = i;
        whereis[i] = i;
    }
    reset();
}

AdvVector::~AdvVector() {
    delete[] data;
    delete[] whereis;
}

void AdvVector::reset() {
    last_id = size;
}

bool AdvVector::thereis( int i ) {
    return whereis[i] < last_id ;
}

bool AdvVector::remove( int i ) {
    i %= size;
    if( whereis[ i ] < last_id ) {
        swap_position( whereis[i], last_id );
        last_id--;
        return true;
    }
    return false;
}

void AdvVector::swap_position( int i, int j ) {
    swap( data[i], data[j] );
    whereis[ data[i] ] = i;
    whereis[ data[j] ] = j;
}

int AdvVector::rnd() {
    if( last_id == 1 )
        return data[ 0 ];
    uniform_int_distribution<int> dist( 0, last_id - 1 );
    return data[ dist( engine ) ];
}

int AdvVector::available_sites() {
    return last_id;
}

bool AdvVector::empty() {
    return last_id == 0;
}