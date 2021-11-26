#include <cstdlib>
#include <iostream>
#include "Grid.h"
#include <fstream>

using namespace std;

int main() {
    Grid g(50);
    g.populate_iid();

    g.print_corr("../data/cf.txt");

    g.print_data("../data/real.txt");

    g.fourier_transform();

    g.print_data("../data/fft.txt", true);

    g.multiply_fft( 0.001 );

    g.print_data("../data/fft2.txt", true);

    g.fourier_transform(true);

    g.print_data("../data/real2.txt");
    
    g.print_corr("../data/cf2.txt");

    system("PAUSE");
}