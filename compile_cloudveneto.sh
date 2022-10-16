rm -f NumericS/main && g++-11 *.cpp -o NumericS/main -std=c++2a -lfftw3 -lm -lgsl -lgslcblas -pthread -march=native -flto -frename-registers -funroll-loops
