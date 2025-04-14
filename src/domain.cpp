#include "domain.h"

//instances of class Rnd.
Rnd rnd;

//domain constructor
Domain::Domain(double x0, double y0, double dx, double dy, int nx, int ny ):x0(x0), y0(y0), dx(dx), dy(dy), nx(nx), ny(ny)
{
    phi = Matrix<double>(nx,ny);
    rho = Matrix<double>(nx,ny);
    efx  = Matrix<double>(nx,ny);
    efy  = Matrix<double>(nx,ny);
    nodes = Matrix<double>(nx,ny);
    Lx = (nx-1)*dx;
    Ly = (ny-1)*dy;
}

void Domain::display(vector<Species> &species)
{
    display::print("\n*********************************************************************************************");
    display::print("         2D Electrostatic Particle-In-Cell Code (ePIC++) with MCC Collision Model ");
    display::print("*********************************************************************************************\n");
    
    display::print(" Simulation Parameters");
    display::print("--------------------------------------------");
    display::print("  System Dimension:           ", nx, " x ", ny);
    display::print("  Area of a Cell:              ",dx*dy);
    display::print("  Normalized Time Step:       ", DT);
    display::print("  Total Time Steps:           ", NUM_TS);
    display::print("  Plasma Density:             ", density);
    
    display::print("\n Plasma Characteristics");
    display::print("--------------------------------------------");
    display::print("  Electron Debye Length:      ", LDe);
    display::print("  Ion Debye Length:           ", LDi);
    display::print("  Electron Plasma Frequency:  ", wpe);
    display::print("  Ion Plasma Frequency:       ", wpi);
    display::print("  Electron Thermal Velocity:  ", LDe * wpe, " m/s");
    display::print("  Ion Thermal Velocity:       ", LDi * wpi, " m/s");
    display::print("  Ion Acoustic Velocity:      ", IAW_vel, " m/s");
    
    display::print("\n Normalization");
    display::print("--------------------------------------------");
    display::print("  Normalization Scheme:       ", 
                        normscheme == 1 ? "Electron-scale" :
                        normscheme == 2 ? "Ion-scale" :
                        normscheme == 3 ? "Sub-cycling" :
                        normscheme == 4 ? "Time in electron, space in ion scale" :
                        "Custom");
    
    display::print("  Velocity Normalization:     ",
                        vel_normscheme == 1 ? "Electron thermal velocity" :
                        vel_normscheme == 2 ? "Ion thermal velocity" :
                        "Ion acoustic velocity");
    
    display::print("  Velocity Normalization Factor:  ", vel_norm);
    display::print("  Simulation Time Period:         ", (2 * Const::PI) / W);
    display::print("  Actual Time Step:               ", DT / W);
    display::print("  Total Simulation Time:          ", (DT / W) * NUM_TS, " seconds");
    display::print("  Ion/Electron Thermal Velocity Ratio: ", vel_ratio);
    display::print("  Electron/Ion Thermal Velocity Ratio: ", 1 / vel_ratio);
    
    display::print("\n Solver and Boundary Conditions");
    display::print("--------------------------------------------");
    display::print("  Boundary Condition:  ", bc);
    display::print("  Potential Solver:    ", SolverType);
    display::print("  Solver Tolerance:    ", tolerance);
    display::print("  Max Iterations:      ", max_iteration);
    display::print("  Shape Function:      ", shapefunction);
    
    display::print("\n Diagnostics");
    display::print("--------------------------------------------");
    display::print("  Diagnostics Type:         ", diagtype);
    display::print("  Write Interval:           ", write_interval);
    //display::print("  Write Interval Phase:     ", write_interval_phase);
    display::print("  Sub-cycle Interval:       ", sub_cycle_interval);
   
    
    display::print("\n Collision Properties");
    display::print("--------------------------------------------");
    display::print("  Neutral Gas Density: ", GAS_DENSITY);
    display::print("  Elastic Collision:   ", enable_elastic_collision ? "✅ Enabled" : "❌ Disabled");
    display::print("  Excitation Collision:", enable_excitation_collision ? "✅ Enabled" : "❌ Disabled");
    display::print("  Ionization Collision:", enable_ionization_collision ? "✅ Enabled" : "❌ Disabled");
    //display::print("  maximum electron collision frequency:", max_electron_coll_freq);
    //display::print("  \u03BD * DT :",max_electron_coll_freq*(DT/wpe));
    
    display::print("\n Execution Mode");
    display::print("--------------------------------------------");
    display::print("  Mode: ", num_threads == 1 ? "🔹 Serial" : "🔹 Parallel (" + std::to_string(num_threads) + " threads)");
    
    int index = 1;
    for(Species &p : species)
    {
        cout << "\n Species (" << index << ") Information ";
        cout << "\n--------------------------------------------";
        cout << "\n  Name:                   " << p.name;
        cout << "\n  Mass:                   " << p.mass;
        cout << "\n  Charge:                 " << p.charge;
        cout << "\n  Temperature:            " << p.temp;
        cout << "\n  Superparticle Weight:   " << p.spwt;
        cout << "\n  Particle Count:         " << p.numparticle;
        cout << "\n  v0x:                     " << p.vsx;
        cout << "\n  v0y:                     " << p.vsy;
        cout << "\n  Normalized Density:     " << p.fract_den;
        cout << "\n  Initialization Type:    " << p.initialization;
        cout << "\n--------------------------------------------\n";
        index++;
    }
    display::print("\n Simulation Ready to Start!\n");
} 

//set normalized parameter.
void Domain:: set_normparam(double LDe, double wpe, double LDi, double wpi)
{
    this->LDe = LDe;
    this->wpe = wpe;
    this->LDi = LDi;
    this->wpi = wpi;  
}

//new code delete old function if this work
void Domain::set_normscheme()
{
    if (normscheme == 1 || normscheme == 3)
    {
        L = LDe;
        W = wpe;
        if(vel_normscheme == 1)
        {
            vel_norm = LDe*wpe;
        }
        if(vel_normscheme == 2)
        {
            vel_norm = LDi*wpi;
        }
        if(vel_normscheme == 3)
        {
            vel_norm = IAW_vel;
        }
        else
        {
            vel_norm = L*W;
        }
        
    }
    else if (normscheme == 2 || normscheme == 4)
    {
        L = LDi;
        W = (normscheme == 2) ? wpi : wpe;  // W remains unchanged for normscheme == 2
        if(vel_normscheme == 1)
        {
            vel_norm = LDe*wpe;
        }
        if(vel_normscheme == 2)
        {
            vel_norm = LDi*wpi;
        }
        if(vel_normscheme == 3)
        {
            vel_norm = IAW_vel;
        }
        else
        {
            vel_norm = L*W;
        }

    }
    else if (normscheme == 5)
    {
        LDe = lenght_scale;
        wpe = time_scale;
        L = LDe;
        W = wpe;
        
        vel_norm = L*W;
    
    }
    else
    {
        throw std::invalid_argument("Invalid normscheme flag");
    }
}

//set time.
void Domain::set_time(double DT, int NUM_TS, int write_interval, int write_interval_phase)
{
    this->DT = DT;
    this->NUM_TS = NUM_TS;
    this->write_interval = write_interval;
    this->write_interval_phase = write_interval_phase;
}

void Domain::ComputeRho(vector<Species> &species)
{
    rho = 0;
    for(auto &sp : species)
    {
        // using overloaded "+" operator for natrix addition and "*" for scalar to matrix multiplication.
        rho += sp.den*(sp.charge*(1/Const::QE));
    }   
}

//scatter densities at grid points
void Domain::Scatter(double lx, double ly,  double value, Matrix<double> &field)
{
    int i = (int)lx;
    int j = (int)ly;
    //cout<<ni<<endl;
    double di = lx-i;
    double dj = ly-j;

    int ngpx = round(lx);
    int ngpy = round(ly);

    if(shapefunction == "NGP")
    {
        field(ngpx,ngpy) += value;
    }

    if(shapefunction == "CIC")
    {
        field(i,j)     += value*(1-di)*(1-dj);
        field(i+1,j)   += value*(di)*(1-dj);
        field(i,j+1)   += value*(1-di)*(dj);
        field(i+1,j+1) += value*(di)*(dj);
    } 
}

//Interpolate density and eletric field value at grid points
double Domain::Interpolate(double lx, double ly, Matrix<double> &field)
{
    int i = (int)lx;
    int j = (int)ly;
    double di = lx - i;
    double dj = ly - j;

    int ngpx = round(lx);
    int ngpy = round(ly);

    double val = 0;

    if(shapefunction == "NGP")
    {
        val = field(ngpx,ngpy);
    }

    if(shapefunction == "CIC")
    {
        val = field(i,j)*(1-di)*(1-dj) + field(i+1,j)*(di)*(1-dj) +
         field(i,j+1)*(1-di)*(dj) + field(i+1,j+1)*(di)*(dj);
    }
    
    return val;
}

double Domain :: XtoL(double pos)
{
    return (pos - x0)/dx;
}
double Domain :: YtoL(double pos)
{
    return (pos - y0)/dy;
}

double Domain::Compute_PE(Species &species)
{
	double pe = 0;
    double unnorm = (density * Const::QE * L) / Const::EPS_0;
	for (int i=0;i < nx;i++)
	{
        for (int j=0;j<ny;j++)
        {
            double ef2 = efx(i,j)*efx(i,j) + efy(i,j)*efy(i,j);
			pe += ef2*(L*L);
        }		
	}
	double pe_norm = 0.5*Const::EPS_0*pe*unnorm*unnorm;

    double Th = (species.temp * Const::eV) * species.spwt * species.numparticle;

    pe_norm /= Th;

    return pe_norm;

}

double Domain::unirand(double lower_bound, double upper_bound) 
{
    // Create a random number generator
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd()); // Seed the generator

    // Define the distribution range
    std::uniform_real_distribution<> distr(lower_bound, upper_bound);

    // Generate and return a random number in the defined range
    return distr(gen);
}
 