#ifndef ADVVECTOR_H
#define ADVVECTOR_H

#include <vector>
#include <random>

class AdvVector {
private:
    int* data;
    int* whereis;
    int last_id;
    static std::minstd_rand engine;

public:
    AdvVector( int size );
    AdvVector(const AdvVector&) = delete; // Copy costruction forbidden
    ~AdvVector();
    void reset();
    int size;

    bool isfree( int i );
    bool remove( int i ); // Return false if it was already removed
    
    int rnd();
    int available_sites();
    bool empty();
};

#endif