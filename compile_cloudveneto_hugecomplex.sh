rm -f CloudVenetoAnalysisHugeButComplex/main && g++-11 *.cpp -o CloudVenetoAnalysisHugeButComplex/main -std=c++2a -lfftw3 -lm -lgsl -lgslcblas -pthread -march=native -flto -frename-registers -funroll-loops
