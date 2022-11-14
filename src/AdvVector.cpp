#include "AdvVector.h"

using namespace std;

minstd_rand AdvVector::engine = minstd_rand{ std::random_device{}() };

AdvVector::AdvVector( int size ) : size( size ) {
    static bool initialized = false;
    if( !initialized ) {
        initialized = true;
        engine.seed( time(NULL) + hash<thread::id>{}(this_thread::get_id()) );
    }
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

bool AdvVector::remove( int i ) {
    i %= size;
    if( whereis[ i ] < last_id ) {
        swap_position( whereis[i], --last_id );
        return true;
    }
    return false;
}

void AdvVector::swap_position( int i, int j ) {
    int temp = data[i];
    data[i] = data[j];
    data[j] = temp;
    whereis[ data[i] ] = i;
    whereis[ data[j] ] = j;
}

int AdvVector::rnd() {
    if( last_id == 1 )
        return data[ 0 ];
    uniform_int_distribution<int> dist( 0, last_id - 1 );
    return data[ dist( engine ) ];
}