#ifndef NEWCFHCALCOLATOR_H
#define NEWCFHCALCOLATOR_H

#include <cstdlib>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <climits>
#include <cstring>
#include <string>
#include "Grid.h"
#include "NewCF.h"

#include <iostream>

namespace NewCFH {

    class Model {
    public:
        std::vector< std::pair<int, std::vector< std::pair<int, int> > > > items;
        std::string keyname;
    };

    class Calculator {
    private:
        NewCF::Model* model;
        Grid<double>* grid;
        std::vector< double >* values;
    public:
        Calculator( NewCF::Model* model, Grid<double>* grid );

        std::vector< double >* calculate();

    };

}

#endif