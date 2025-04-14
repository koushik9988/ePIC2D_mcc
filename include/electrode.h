#ifndef ELECTRODE_H
#define ELECTRODE_H

#include "domain.h"
#include <iostream>

class Domain; // Forward declaration of Domain class

class Electrode
{
    public:
    double x0, y0;    // Center of the electrode
    double r;         // Radius
    double voltage;   // Voltage

    Electrode(double x0, double y0, double r, double voltage, Domain& domain);
    bool Iselectrode(double x, double y);
    private:
    Domain &domain;
};

#endif