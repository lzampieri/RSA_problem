#ifndef POLYMER_H
#define POLYMER_H

#include <vector>
#include <random>
#include <string>
#include "Grid.h"

#ifndef INT_MAX
#include <climits>
#endif

class Polymer : public GridProps
{
public:
    std::vector<GridSite> *atoms;
    std::vector<GridSite> *blocked; // Sites in which cannot be placed a polymers, considering that (0,0) is occupied

    Polymer(std::vector<GridSite> gss);
    ~Polymer();

    bool canStay(const Grid<int> &grid, const GridSite position);
};

class Polymers
{
public:
    std::vector<Polymer *> variants;
    int N;
    std::string keyname;

    Polymers(std::string keyname);
    ~Polymers();

    void addVariant(Polymer *p);
    inline Polymer* operator[](int i) { return variants[i]; }
    inline Polymer* at(int i) { return variants[i]; };

    Polymers(const Polymers &) = delete; // Copy costruction forbidden
};

class PolymersFactory
{
public:
    virtual std::string factname() const { return "Error"; }
    virtual Polymers *create(GridProps &gp) = 0;
    Polymers *operator()(GridProps &gp) { return create(gp); }
    Polymers *create(int s)
    {
        GridProps gp(s);
        return create(gp);
    }
    Polymers *operator()(int s) { return create(s); }

    static std::vector<PolymersFactory *> StdPolymers;
};

namespace StdPolymers
{
    class Dimers : public PolymersFactory
    {
        std::string factname() const { return "Dimers"; };
        Polymers *create(GridProps &gp);
    };
    class LinearTrimers : public PolymersFactory
    {
        std::string factname() const { return "LinearTrimers"; }
        Polymers *create(GridProps &gp);
    };
    class Trimers : public PolymersFactory
    {
        std::string factname() const { return "Trimers"; }
        Polymers *create(GridProps &gp);
    };
    class Squared : public PolymersFactory
    {
        std::string factname() const { return "Squared"; }
        Polymers *create(GridProps &gp);
    };
    class LinearPentamers : public PolymersFactory
    {
        std::string factname() const { return "LinearPentamers"; }
        Polymers *create(GridProps &gp);
    };
}

#endif