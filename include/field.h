#ifndef _FIELD_H_
#define _FIELD_H_

#include "domain.h"
#include <cstring>
#include "linalg.h"
#include "function.h"
#include "solvers.h"
#include "grid.h"

class Domain; // Forward declaration of Domain class
class Grid; // Forward declaration of Grid class

class FieldSolve
{
public:
    // Constructor
    //FieldSolve(Domain &domain,Grid *grid):domain(domain),grid(grid){};
    FieldSolve(Domain &domain, const std::vector<Grid*> &grids) : domain(domain), grids(grids) {}

    void nrpcgsolver();
    void pcgsolver();
    bool gssolver();
    bool gsfluidsolver();
    void CalculateEfield();
    void PotentialSolver();
    void spectral();

private:
    Domain &domain;
    std::vector<Grid*> grids; // Vector of pointers to Grid objects
};

#endif 
