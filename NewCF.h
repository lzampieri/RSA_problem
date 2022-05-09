#ifndef NEWCORRFUNCCALCOLATOR_H
#define NEWCORRFUNCCALCOLATOR_H

#include <cstdlib>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <climits>
#include <cstring>
#include <string>
#include "Grid.h"

#include <iostream>

namespace NewCF {

    class Model {
    public:
        std::vector< std::pair<int, std::vector< std::pair<int, int> > > > items;
        std::string keyname;
    };

    class Calculator {
    private:
        Model* model;
        Grid<int>* grid;
        std::vector< double >* values;
    public:
        Calculator( Model* model, Grid<int>* grid );

        std::vector< double >* calculate();

    };

    // class Equispaced : public Model {
    // public:
    //     Equispaced(int start, int end, int step);
    // };

    class Expospaced : public Model {
    public:
        Expospaced(double base, int maxval);
    };
}

#endif