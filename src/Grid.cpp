#include "Grid.h"

using namespace std;

// GridProps
// Constructors
GridProps::GridProps(int d1, int d2) : d1(d1), d2(d2), d1d2(d1 * d2) {}
GridProps::GridProps(int d) : d1(d), d2(d), d1d2(d * d) {}
GridProps::GridProps(const GridProps &G) : d1(G.d1), d2(G.d2), d1d2(G.d1 * G.d2) {}

// Distances

double GridProps::d(int x1, int y1, int x2, int y2) const
{
    double deltaX = _pbc_x(x1 - x2);
    deltaX = deltaX < d1 - deltaX ? deltaX : d1 - deltaX;
    double deltaY = _pbc_y(y1 - y2);
    deltaY = deltaY < d2 - deltaY ? deltaY : d2 - deltaY;
    return sqrt(deltaX * deltaX + deltaY * deltaY);
}

// GridSite
// Constructors

GridSite::GridSite(int X, int Y, const GridProps &gp) : GridProps(gp)
{
    setX(X);
    setY(Y);
}
GridSite::GridSite(int I, const GridProps &gp) : GridSite(gp._xy(I)) {}
GridSite::GridSite(const GridSite &G) : GridProps(G)
{
    setX(G.X);
    setY(G.Y);
}

// Operators
GridSite &GridSite::operator=(const GridSite &gs)
{
    d1 = gs.d1;
    d2 = gs.d2;
    d1d2 = gs.d1 * gs.d2;
    setX(gs.X);
    setY(gs.Y);
    return *this;
}

GridSite &GridSite::operator+=(const GridSite &xy2)
{
    setX(X + xy2.X);
    setY(Y + xy2.Y);
    return *this;
}

GridSite &GridSite::operator-=(const GridSite &xy2)
{
    setX(X - xy2.X);
    setY(Y - xy2.Y);
    return *this;
}

// Grid
// Constructors

template <class T>
Grid<T>::Grid(int d) : Grid<T>(d, d){};

template <class T>
Grid<T>::Grid(int d1, int d2) : GridProps(d1, d2)
{
    u = new vector<T>( d1d2, 0);
}

template <class T>
Grid<T>::~Grid()
{
    u->clear();
    delete u;
}

template <class T>
void Grid<T>::print_data(const string filename) const
{
    ofstream out(filename);
    for (int i = 0; i < d1; i++)
        for (int j = 0; j < d2; j++)
        {
            out << i << '\t' << j << '\t' << operator()(i, j) << '\n';
        }
    out.close();
}