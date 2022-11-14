#ifndef ADVVECTOR_H
#define ADVVECTOR_H

#include <vector>
#include <random>
#include <time.h>
#include <thread>

class AdvVector {
public:
    int* data;
    int* whereis;
    int last_id;
    static std::minstd_rand engine;

public:
    AdvVector( int size );
    AdvVector(const AdvVector&) = delete; // Copy costruction forbidden
    ~AdvVector();
    inline void reset() { last_id = size; }
    int size;

    inline bool thereis( int i ) const { return whereis[i] < last_id; };
    bool remove( int i ); // Return false if it was already removed

    void swap_position( int i, int j );
    
    int rnd();
    inline int available_sites() const { return last_id; };
    inline bool empty() const { return last_id == 0; }
};

#endif