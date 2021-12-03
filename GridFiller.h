#ifndef GRIDFILLER_H
#define GRIDFILLER_H

#include <ctime>

class Grid;

class GridFiller {

public:
    static void iid(Grid& g);
    static void coscos(Grid& g, double fact1, double fact2);
    static void square(Grid& g, int side);
};


#endif