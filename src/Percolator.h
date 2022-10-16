#ifndef PERCOLATOR_H
#define PERCOLATOR_H

#include <vector>
#include "Grid.h"

class Percolator {
    // Not multithreading by choice,
    // considering the multithreading
    // is already implemented in the
    // replicas
private:
    Grid<int>& grid;
    std::vector<bool> inserted;
    std::vector<GridSite> to_check;
    int to_check_end; // The place just after the last occupied one

    bool next(GridSite& i);
    bool add(GridSite i);
    void reset();
public:
    Percolator(Grid<int>& grid);

    bool is_percolating(short what);
};

#endif