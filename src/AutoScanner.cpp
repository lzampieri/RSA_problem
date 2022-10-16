#include "AutoScanner.h"

using namespace std;

string AutoScanner::populate() {
    string folder = computeFolder();

    for( int s : sides )
        for( double g : gammas )
            for( double q : qs )
                for( PolymersFactory* p : ps )
                    rps.push_back( ReplicatorParams(
                     // Size DefectsFracs Gamma ChunkSize   Tolerance  CFModel  Polymers                                   Percolation  NThreads,  Draw  Verbose  SavePath
                        s,   q,           g,    chunk_size, tolerance, CFmodel, ( p == nullptr ? nullptr : p->create(s) ), percolation, n_threads, draw, verbose, folder
                    ));

    return folder;
}

void AutoScanner::loadFromTxt( string file ) {
    ifstream in( file );
    string temp;

    in  >> chunk_size
        >> tolerance;
    
    smart_getline( in, temp );
    string_to_array( temp, sides );
    
    smart_getline( in, temp );
    string_to_array( temp, gammas );

    smart_getline( in, temp );
    string_to_array( temp, qs );

    smart_getline( in, temp );
    string_to_array( temp, ps );
    
    in >> percolation
        >> draw
        >> verbose
        >> n_threads;

    filename = file;
}

string AutoScanner::computeFolder( ) {

    string folder = filename;

    // Ensure only numbers and letters
    folder.erase( remove_if( folder.begin(), folder.end(),
        [](char c) { return !isalpha(c) && !isdigit(c); } ),
        folder.end()
        );

    // Append current date
    folder += "_" + date::format("%Y%m%d", chrono::system_clock::now());

    // Check for existence
    if( filesystem::exists( folder ) ) {
        int progressive = 0;
        // If already exists, append progressive
        while( filesystem::exists( folder + "_" + to_string( progressive ) ) )
            progressive ++;
        folder += "_" + to_string( progressive );
    }

    // Create
    filesystem::create_directory( folder );

    // Save there a copy of the configuration file
    saveToTxt( folder + "/params.scan" );

    return folder;
}

void AutoScanner::saveToTxt( string file ) {
    ofstream out( file );

    out << chunk_size << '\n'
        << tolerance << '\n'
        << array_to_string( sides ) << '\n'
        << array_to_string( gammas ) << '\n'
        << array_to_string( qs ) << '\n'
        << array_to_string( ps ) << '\n'
        << percolation << '\n'
        << draw << '\n'
        << verbose << '\n'
        << n_threads << endl;

    out.close();
}

void AutoScanner::string_to_array( const string text, vector<double>& vec ) {
    stringstream ss( text );
    vec.clear();
    double t;
    while( ss >> t )
        vec.push_back( t );
}

void AutoScanner::string_to_array( const string text, vector<int>& vec ) {
    stringstream ss( text );
    vec.clear();
    int t;
    while( ss >> t )
        vec.push_back( t );
}

void AutoScanner::string_to_array( const string text, vector<PolymersFactory*>& vec ) {
    stringstream ss( text );
    string temp;
    vec.clear();
    while( ss >> temp ) {
        for( PolymersFactory* p : PolymersFactory::StdPolymers ) {
            if( p->factname() == temp ) {
                vec.push_back( p );
                break;
            }
        }
    }
}

string AutoScanner::array_to_string( const std::vector<double>& vec ) {
    stringstream ss( "" );
    for( double t : vec )
        ss << t << ' ';
    return ss.str();
}

string AutoScanner::array_to_string( const std::vector<int>& vec ) {
    stringstream ss( "" );
    for( int t : vec )
        ss << t << ' ';
    return ss.str();
}

string AutoScanner::array_to_string( const std::vector<PolymersFactory*>& vec ) {
    stringstream ss( "" );
    for( const PolymersFactory* t : vec )
        ss << ( t == nullptr ? "Null ptr" : t->factname() ) << ' ';
    return ss.str();
}

void AutoScanner::smart_getline( ifstream& in, string& str ) {
    do {
        getline( in, str );
    } while( str.length() == 0 );
}