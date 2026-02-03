#ifndef _SPECIES_H_
#define _SPECIES_H_

#include <iostream>
#include <list>
#include <vector>
#include "domain.h"
#include <thread>
#include <cstring>
#include <cmath>
#include <iostream>
#include <vector>
#include "slap.h"
#include <functional>
#include "grid.h"

using namespace std;

class Domain;
class Grid;

class Particle
{
    public:
    double pos[2];
    double vel[3];

    //particle constructor
    Particle(double x, double y, double vx, double vy, double vz)
    {
        pos[0] = x;
        pos[1] = y;
        vel[0] = vx;
        vel[1] = vy;
        vel[2] = vz;
    }
};

class Species
{
    public:
    /// species name
    string name;
    /// species mass
    double mass;
    /// species charge
    double charge;
    /// @brief specific weight
    double spwt;
    /// @brief species temparature
    double temp;
    /// @brief density
    Matrix<double> den;
    /// @brief no of simulation particle
    int numparticle;
    /// @brief 
    int charge_sig;
    std::string initialization_posx;
    std::string initialization_posy;
    std::string initialization_velx;
    std::string initialization_vely;
    double vsx;
    double vsy;
    double fract_den;

    vector<Particle> part_list;

    Species(string name, double mass, double charge, double spwt, double temp, int numparticle, double vsx, double vsy, double fract_den , 
        std:: string initialization_posx, std:: string initialization_posy,std:: string initialization_velx, std:: string initialization_vely, Domain &domain);

    //declare member functions or methods
    void AddParticle(Particle part);
    void Push_species_serial(int sub_cycle, const std::vector<Grid*> &grids);
    void ScatterSpecies_serial();
    void Rewind_species();
    vec<double> Compute_KE(Species &species);
    vec<double> Compute_Momentum(Species &species);
    
    private:
    Domain &domain;  
};



//new
struct SpeciesParams
{
    string name;
    double mass;
    int num;
    double temp;
    int charge_sign;
    double normden;
    double vsx;
    double vsy;
    string loadtype_posx;
    string loadtype_posy;
    string loadtype_velx;
    string loadtype_vely;
};


#endif 