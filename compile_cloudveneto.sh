rm -f CloudVenetoAnalysis/main && g++-11 *.cpp -o CloudVenetoAnalysis/main -std=c++2a -lfftw3 -pthread -march=native -flto -frename-registers -funroll-loops
