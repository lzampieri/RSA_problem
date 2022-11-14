rm -f NewDepositionData/main && g++-11 src/*.cpp -o NewDepositionData/main -std=c++2a -lfftw3 -lm -lgsl -lgslcblas -pthread -march=native -flto -frename-registers -funroll-loops
