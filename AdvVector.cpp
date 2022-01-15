#include "AdvVector.h"

using namespace std;

minstd_rand AdvVector::engine = minstd_rand{ std::random_device{}() };

AdvVector::AdvVector( int size ) : size( size ) {
    data = new int[size];
    whereis = new int[size];
    reset();
}

void AdvVector::reset() {
    for( int i=0; i < size; i++ ) {
        data[i] = i;
        whereis[i] = i;
    }
    last_id = size;
}

bool AdvVector::isfree( int i ) {
    return whereis[i] >= 0;
}

bool AdvVector::remove( int i ) {
    if( whereis[i] >= 0 ) {
        data[ whereis[i] ] = data[ --last_id ];
        whereis[ data[ whereis[i] ] ] = whereis[i];
        whereis[i] = -1;
        return true;
    }
    return false;
}

int AdvVector::rnd() {
    uniform_int_distribution<int> dist( 0, last_id - 1 );
    return data[ dist( engine ) ];
}

int AdvVector::available_sites() {
    return last_id;
}

bool AdvVector::empty() {
    return last_id == 0;
}