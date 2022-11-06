#ifndef AUTOSCANNER_H
#define AUTOSCANNER_H

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <filesystem>

#include "Replicator.h"
#include "Polymer.h"
#include "NewCF.h"
#include "date.h"

class AutoScanner {
    
private:

    int chunk_size;
    double tolerance;

    std::vector<int> sides;
    std::vector<double> gammas;
    std::vector<double> qs;
    std::vector<PolymersFactory*> ps = { nullptr };

    bool draw;

public:
    bool verbose;

private:
    int n_threads;

    std::string filename;

    NewCF::Model* CFmodel;

public:
    std::vector< ReplicatorParams > rps;

public:
    AutoScanner();

    std::string populate();

    void loadFromTxt( std::string file );
    void saveToTxt( std::string file );

private:
    std::string computeFolder( );

    void string_to_array( const std::string text, std::vector<double>& vec );
    void string_to_array( const std::string text, std::vector<int>& vec );
    void string_to_array( const std::string text, std::vector<PolymersFactory*>& vec );
    
    std::string array_to_string( const std::vector<double>& vec );
    std::string array_to_string( const std::vector<int>& vec );
    std::string array_to_string( const std::vector<PolymersFactory*>& vec );

    bool smart_getline( std::ifstream& in, std::string& str );
};

#endif