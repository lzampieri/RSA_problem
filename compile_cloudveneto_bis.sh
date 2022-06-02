rm -f CloudVenetoAnalysisBis/main && g++-11 *.cpp -o CloudVenetoAnalysisBis/main -std=c++2a -lfftw3 -lm -lgsl -lgslcblas -pthread -march=native -flto -frename-registers -funroll-loops
