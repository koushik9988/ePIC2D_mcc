// header file(.h or .hpp) contain function,class,namespace etc declaration
#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "species.h"
#include <iostream>
#include <vector>
#include <random>
#include <stdexcept>
#include "output.h"
#include <cstring>
#include "math.h"
#include "slap.h"


class Species;

using namespace std;

/**
 * @class Domain
 * @brief Represents the simulation domain for the 1D Electrostatic PIC (Particle-in-Cell) simulation.
 */

/// @brief 
class Domain
{
    public:
    ///origin points
    double x0; 
    double y0;
    ///max lenght of system
    double Lx;
    double Ly;
    ///cell spacing 
    double dx;
    double dy;
    ///no of cell 
    int cell_no;
    ///no of grid points
    int nx;
    int ny;
    ///total simulation time step 
    int NUM_TS;
    ///interval for writing field data
    int write_interval;
    ///interval for writing phase space data
    int write_interval_phase;
    ///interval for sub_cycling ion motion
    int sub_cycle_interval;
    //double LD;

    ///data structures to hold data
    Matrix<double> phi;
    Matrix<double> rho;
    Matrix<double> efx;
    Matrix<double> efy;

    Matrix<double> nodes;

    ///electron debye lenght
    double LDe;
    ///electron plasma frequency
    double wpe;
    ///Ion debye lenght
    double LDi;
    ///ion plasma frequency
    double wpi;
    ///user defined lenght scale
    double lenght_scale;
    ///user defined time scale
    double time_scale;
    ///normalization scale
    double W,L;
    ///Ion acoustic speed
    double IAW_vel;
    ///velocity normalization factor
    double vel_norm;

    double density;
    ///time step coefficient
    double DT;
    

    ///no of species
    int species_no;
    ///no of electron crossed more than one cell
    int ele_cross = 0;
    ///no of ion crossed more than one cell
    int ion_cross = 0;
    ///average no of cell crossed by eelctron
    int crossed_cellno_ele = 0;
    ///average no of cell crossed ny ion
    int crossed_cellno_ion = 0;
    ///flag for normalization scheme
    double norm;
    /// name of potential solver.
    std::string SolverType ;
    /// tolarance of iterative solver
    double tolerance;
    /// @brief  maximum no of iteration for interative solver
    int max_iteration;
    ///normalizing velocty
    //double vel_norm;//= domain.L*domain.W;
    double vel_ratio;
    ///flag to on/off parrticle pusher parallization
    bool push_parallal;
    ///flag to on/off charge deposition parallization
    bool deposit_parallal;
    ///number of cpu threads to run the simulation.
    int num_threads;

    double error;

    std::string pload;
    
    std::string bc;

    std::string shapefunction;

    std::string diagtype;

    int normscheme;
    int vel_normscheme;

    double delta_g;

    //new
    bool enable_elastic_collision;    // Flag for elastic collisions
    bool enable_excitation_collision; // Flag for excitation collisions
    bool enable_ionization_collision; // Flag for ionization collisions
    double GAS_DENSITY;


    /**
     * @brief Constructor for the Domain class.
     * @param x0 Origin point.
     * @param dx Cell spacing.
     * @param ni Number of grid points.
     */
    Domain(double x0, double y0, double dx, double dy, int nx, int ny);
 

    void display(vector<Species> &species);

    /**
     * @brief function to set normalization paramter.
     * @param LD debye lenght.
     * @param wp plasma frequency.
     * @param W normalizing frequency.
     */
    void set_normparam(double LDe, double wpe, double LDi, double wpi);
    /**
     * @brief function to set normalization scheme.
     * @param normscheme string argument to decide normalization scheme.
     */
    void set_normscheme();
    /**
     * @brief function to set different simulation time parameter.
     * @param DT time step.
     * @param NUM_TS Total simulation time step.
     * @param write_interval data writing interval.
     */
    void set_time(double DT, int NUM_TS, int write_interval, int write_interval_phase);
    /**
     * @brief function to Compute Charge Density.
     * @param species Species class instance.
     */
    void ComputeRho(vector<Species> &species);
    /**
     * @brief function to scatter density to mesh/grid.
     * @param lc co-ordinate of particle
     * @param value specific weight of particle
     * @param field species density
     */ 
    void Scatter(double lx, double ly,  double value, Matrix<double> &field);
    /**
     * @brief function to calculate field value for each particle depending on their location.
     * @param lc co-ordinate of particle
     * @param field field value(electric field)
     */
    double Interpolate(double lx, double ly, Matrix<double> &field);
    /**
     * @brief function convert physical coordinate into Logical coordinate.
     * @param pos physical loaction of the charged particle.
     */
    double XtoL(double pos);
    double YtoL(double pos);
    /**
     * @brief function to calculate electric field potential energy.
     * @param species name of the normalizing species.
     */
    /// @brief structures to hold density value for each cpu threads .
    std::vector<vec<double>> buffers;

    double Compute_PE(Species &species);

    double unirand(double lower_bound, double upper_bound);

    int ionfixed;
};

//define a namespace name Const
namespace Const
{
	const double EPS_0 = 8.85418782e-12;  	//vacuum permittivity
	const double QE = 1.602176565e-19;		//electron charge
	const double AMU = 1.660538921e-27;		//atomic mass unit
	const double ME = 9.10938215e-31;		//electron mass
	const double K_b = 1.380648e-23;			//Boltzmann constant
	const double PI = 3.141592653;			//pi
	const double EV_to_K = QE/K_b;				//1eV in K ~ 11604
    const double eV = 1.602176565e-19;		//1 eletron volt
    
}

/*object for sampling random numbers*/
class Rnd
{
    public:
	//constructor: set initial random seed and distribution limits
	//Rnd(): mt_gen{std::random_device()()}, rnd_dist{0,1.0} {}
    Rnd(): mt_gen{0}, rnd_dist{0,1.0} {}
	double operator() () {return rnd_dist(mt_gen);}

protected:
	std::mt19937 mt_gen;	    //random number generator
	std::uniform_real_distribution<double> rnd_dist;  //uniform distribution
};

extern Rnd rnd;		//tell the compiler that an object of type Rnd called rnd is defined somewhere

#endif