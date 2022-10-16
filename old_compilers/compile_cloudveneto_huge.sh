rm -f CloudVenetoAnalysisHuge/main && g++-11 *.cpp -o CloudVenetoAnalysisHuge/main -std=c++2a -lfftw3 -lm -lgsl -lgslcblas -pthread -march=native -flto -frename-registers -funroll-loops
