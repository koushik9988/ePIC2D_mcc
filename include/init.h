#ifndef _INIT_H_
#define _INIT_H_

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "species.h"
#include "domain.h"
#include "output.h"
#include "slap.h"
#include "iniparser.h"

class Domain;
class Species;

class Init 
{
public:
    //in the constructor we pass Species and Domain class instances to access the species and domain class
    //memebr variable and methods 
    Init(Species &species, Domain &domain);
    //void inject(Species &species, int numParticles);
    static double SampleVel(Species &species);
    static double SampleVel(Species &species, double temp);
    
    //void inject(Species &species, Domain &domain, int numParticles);

private:
    Species &species;
    Domain &domain;
    //double* SampleVel(Species &species);
};

void inject(Species &species, Domain &domain, int numParticles);

#endif 
