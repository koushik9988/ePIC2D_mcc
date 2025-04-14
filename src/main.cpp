/* */
#include "iniparser.h"
#include <iostream>
#include "field.h"
#include "init.h"
#include "domain.h"
#include "species.h"
#include <cmath>
#include <fstream>
#include <time.h>
#include <chrono>
#include "output.h"
#include "extrafun.h"
#include <thread>
#include "slap.h"
#include "collision.h"
#include "electrode.h"
#include "grid.h"
#include "emitter.h"

using namespace std; 
using namespace display;

int main( int argc , char *argv[]) 
{     
    if(argc<2)
    {
      cout<<"ERROR, at least one argument expected (the input file)."<<endl;
      exit (EXIT_FAILURE);
    }

    print("---------------------------------1D Electrostatic Particle In Cell code-------------------------------------");
    //parsing input.ini file and storing values
    const std::string filename = argv[1];
    
    auto iniData = INIParser::parse(filename);

    auto species_section = iniData["species"];
    int species_no = species_section.size();

    //file output
    std::string outputfolder = INIParser::getString(iniData["file"],"output");
    //grid or domain
    int Nx = INIParser::getInt(iniData["domain"], "Nx");//cell no
    int Ny = INIParser::getInt(iniData["domain"], "Ny");//cell no
    int nx = Nx+1; // no of grid points is one more than cell no
    int ny = Ny+1;
    double x0 = INIParser::getDouble(iniData["domain"],"x0");
    double y0 = INIParser::getDouble(iniData["domain"],"y0");

    
    //collision
    std::string elastic_flag = INIParser::getString(iniData["collision"],"elastic");
    std::string excitation_flag = INIParser::getString(iniData["collision"],"excitation");
    std::string ionization_flag = INIParser::getString(iniData["collision"],"ionization");
    double GAS_DENSITY = INIParser::getDouble(iniData["collision"],"GAS_DENSITY");
   
    //time
    int NUM_TS = INIParser::getInt(iniData["time"], "NUM_TS");
    double DT_coeff = INIParser::getDouble(iniData["time"],"DT_coeff");

    //simulation parameter
    std::string shapefunction = INIParser::getString(iniData["simulation"],"shapefunction");
    double den = INIParser::getDouble(iniData["simulation"],"density");
    std:: string bc = INIParser::getString(iniData["simulation"],"bc");
    std::string push_parallal = INIParser::getString(iniData["simulation"],"push_parallal");
    std::string deposit_parallal = INIParser::getString(iniData["simulation"],"deposit_parallal");
    int ionfixed = INIParser::getInt(iniData["simulation"],"ionfixed");
    

    //Diagnostics
    int precision = INIParser::getInt(iniData["diagnostics"],"precision");
    int write_flag = INIParser::getInt(iniData["diagnostics"],"write_flag");
    int write_interval = INIParser::getInt(iniData["diagnostics"],"write_interval");
    int write_interval_phase = INIParser::getInt(iniData["diagnostics"],"write_interval_phase");
    int write_diagnostics = INIParser::getInt(iniData["diagnostics"],"write_diagnostics");
    int save_fig = INIParser::getInt(iniData["diagnostics"],"save_fig");
    int sub_cycle_interval = INIParser::getInt(iniData["diagnostics"],"sub_cycle_interval");
    

    //normalization
    int norm_scheme = INIParser::getInt(iniData["normalization"],"norm_scheme");
    int vel_normscheme = INIParser::getInt(iniData["normalization"],"vel_norm_scheme");
    double L_scale = INIParser::getDouble(iniData["normalization"],"lenght_scale");
    std::string T_scale = INIParser::getString(iniData["normalization"],"time_scale");


    std::string SolverType = INIParser::getString(iniData["solver"],"solvertype");
    double tolerance = INIParser::getDouble(iniData["solver"],"tolerance");
    double max_iteration = INIParser::getInt(iniData["solver"],"max_iteration");

    
    //vector to store species data
    std::vector<std::string> names;
    std::vector<double> mass;
    std::vector<int> nParticles;
    std::vector<double> temps;
    std::vector<int> charge_signs;
    std::vector<double> frac_densities;
    std::vector<double> normden;
    std::vector<double> vsx;
    std::vector<double> vsy;
    std::vector<std::string> pos_init;

    print(species_no);

    names.reserve(species_no);
    mass.reserve(species_no);
    nParticles.reserve(species_no);
    temps.reserve(species_no);
    charge_signs.reserve(species_no);
    frac_densities.reserve(species_no);
    normden.reserve(species_no);
    vsx.reserve(species_no);
    vsy.reserve(species_no);
    pos_init.reserve(species_no);

    // First pass: Parse the species section and populate vectors
    for (const auto& species_entry : species_section)
    {
        const std::string& line = species_entry.second;
        std::vector<std::string> tokens = INIParser::split(line, ',');

        if (tokens.size() == 9)
        {  
            names.push_back(tokens[0]);                      // Species name
            mass.push_back(std::stod(tokens[1]));
            nParticles.push_back(std::stoi(tokens[2]));      // Number of particles
            temps.push_back(std::stod(tokens[3]));           // Temperature
            charge_signs.push_back(std::stoi(tokens[4]));    // Charge sign (integer -1 or 1)
            frac_densities.push_back(std::stod(tokens[5]));  // Fractional density
            vsx.push_back(std::stod(tokens[6])); 
            vsy.push_back(std::stod(tokens[7])); 
            pos_init.push_back(tokens[8]);
        }
    }

    double k = 0;
    for(int i = 0 ; i < species_no; i++)
    {
        k += (-charge_signs[i])*frac_densities[i];
    }

    //display::print(k);

    normden[1] = den;
    normden[0] = den/k;
    
    for(int i = 2 ;i < species_no; i++)
    {
        normden[i] = frac_densities[i]*normden[0];
    }

    //noramlizing quantity(electron)
    double LDe = sqrt((Const::EPS_0*Const::K_b*temps[0]*Const::EV_to_K)/(normden[0]*Const::QE*Const::QE)); // Electron Debye Length   
    double wpe = sqrt((normden[0]*Const::QE*Const::QE)/(mass[0]*Const::EPS_0)); // Total Electron Plasma Frequency
    double wpi = sqrt((normden[1]*Const::QE*Const::QE)/(mass[1]*Const::EPS_0)); //ion timescale
    double LDi = sqrt((Const::EPS_0*Const::K_b*temps[1]*Const::EV_to_K)/(normden[1]*Const::QE*Const::QE)); // ion Debye Length
    double CS = sqrt(temps[0]*Const::K_b*Const::EV_to_K/mass[1]); // Ion acoustic speed

    //print("electron plasma freq",wpe);

    double vthe = LDe*wpe;
    double vthi = LDi*wpi;

    //user defined scales
    double energy_scale = 1;
    double time_scale;
    double lenght_scale = L_scale;
    if(T_scale == "omegape")
    {
        time_scale = wpe; 
    }
    else if(T_scale == "omegapi")
    {
        time_scale = wpi; 
    }
    else if(T_scale != "omegape" || T_scale != "omegapi")
    {
        time_scale = std::stod(T_scale);
    }
 
    double dx;
    double dy;
    double DT;
    double stepSize;
    if(norm_scheme == 2 || norm_scheme == 4)
    {
        stepSize = LDi;
        dx = stepSize/LDi;
        dy = stepSize/LDi;
        DT = DT_coeff*(1.0/wpi);
        DT = wpi*DT;
    }
    if(norm_scheme == 1 || norm_scheme == 3)
    {
        stepSize = LDe;
        dx = stepSize/LDe;
        dy = stepSize/LDe;
        DT = DT_coeff*(1.0/wpe);
        DT = wpe*DT;
    }
    
    if(norm_scheme == 5)
    {
        stepSize = lenght_scale;
        dx = stepSize/lenght_scale;
        dy = stepSize/lenght_scale;
        DT = DT_coeff*(1.0/time_scale);
        DT = time_scale*DT;
    }
    
    Domain domain(x0,x0,dx,dy,nx,ny);
    
    domain.bc  = bc;
    domain.SolverType = SolverType;
    domain.tolerance = tolerance;
    domain.species_no = species_no;
    domain.density = den;
    domain.shapefunction = shapefunction;
    domain.IAW_vel = CS;
    //domain.num_threads = num_threads;
    //domain.push_parallal = string_to_bool(push_parallal);
    //domain.deposit_parallal = string_to_bool(deposit_parallal);
    //domain.pload = pload;  
    domain.max_iteration = max_iteration;
    domain.normscheme = norm_scheme;
    domain.sub_cycle_interval = sub_cycle_interval;
    domain.set_normparam(LDe,wpe,LDi,wpi);
    domain.set_time(DT,NUM_TS,write_interval,write_interval_phase);
    domain.set_normscheme();
    domain.ionfixed = ionfixed;
    domain.enable_elastic_collision = string_to_bool(elastic_flag);    // Flag for elastic collisions
    domain.enable_excitation_collision = string_to_bool(excitation_flag); // Flag for excitation collisions
    domain.enable_ionization_collision = string_to_bool(ionization_flag); // Flag for ionization collisions
    domain.GAS_DENSITY = GAS_DENSITY;
    


    std::vector<Species> species_list;
    species_list.reserve(species_no);
 
    for (int i = 0 ;i < species_no; i++)
    {
        double computed_spwt = 0;

        computed_spwt = normden[i] * domain.Lx * domain.Ly * domain.L *domain.L / nParticles[i];
        //print(domain.L);

        species_list.emplace_back(names[i], mass[i], charge_signs[i]*Const::QE, computed_spwt, temps[i], 
            nParticles[i], vsx[i], vsy[i], frac_densities[i], pos_init[i], domain);
    }

    domain.display(species_list);

    print("Press Enter to continue...");
    std::cin.get();

    Output output(outputfolder,domain);

    //output.write_metadata(Nx, Ny, NUM_TS, write_interval, write_interval_phase, DT, nE, nI, tempE, tempI, 
    //massI,  den, save_fig, vxe, vye, vxi, vyi, domain.normscheme, domain.sub_cycle_interval);

    output.precision = precision;
    //output.Energy_plot =  Energyplot_flag ;
    //output.Potentialfield_plot = Potentialfieldplot_flag;
    //output.Chargedensity_plot = chargeplot_flag ;
    //output.keflag = keflag;
    //output.peflag = peflag;
    //output.teflag = teflag;
    //output.phase_plot = phaseplot_flag;
    //output.dft_flag = dft_flag;
    //output.species_index = species_index;
    output.write_metadata(nx,NUM_TS,write_interval,write_interval_phase,DT_coeff,den,save_fig,domain.normscheme,
        domain.sub_cycle_interval,LDe,LDi,wpe,wpi,species_no,0);
    output.write_species_metadata(species_list);


    //output.write_metadata(NC,NUM_TS,write_interval,write_interval_phase,DT_coeff, nE, nI, nN, nB,tempE,tempI,tempB, tempN, alpha, beta, massI, massN, massB,den,save_fig,v_e, v_i,v_n,v_b,domain.normscheme,domain.sub_cycle_interval);
    
    auto start_time = std::chrono::high_resolution_clock::now();

    ////-----grid initialization------------------
    //new code 

    std::vector<Grid*> grids;
    const auto& gridsSection = iniData.at("Grids");

    // Get the number of grids
    int grid_count = INIParser::getInt(gridsSection, "grid_count");

    // Loop through the number of grids
    for (int i = 0; i < grid_count; i++)
    {
        std::string prefix = "grid_" + std::to_string(i) + ".";

        std::string type = INIParser::getString(gridsSection, prefix + "type");

        GridParams grid_params;
        RectangularGridParams grid_params1;

        if (type == "circular")
        {
            grid_params.grid_centerx    = INIParser::getInt(gridsSection, prefix + "x_center");
            grid_params.grid_centery    = INIParser::getInt(gridsSection, prefix + "y_center");
            grid_params.num_electrodes  = INIParser::getInt(gridsSection, prefix + "electrode_number");
            grid_params.grid_rad        = INIParser::getDouble(gridsSection, prefix + "grid_radius");
            grid_params.electrode_rad   = INIParser::getDouble(gridsSection, prefix + "electrod_radius");
            grid_params.voltage         = (Const::eV / (Const::K_b * Const::EV_to_K)) * INIParser::getDouble(gridsSection, prefix + "electrode_voltage");

            grids.push_back(Create_Grid(GridType::Circular, domain, grid_params, grid_params1));
        }
    
        else if (type == "reactconduct")
        {
            grid_params1.min_x = INIParser::getInt(gridsSection, prefix + "min_x");
            grid_params1.min_y = INIParser::getInt(gridsSection, prefix + "min_y");
            grid_params1.max_x = INIParser::getInt(gridsSection, prefix + "max_x");
            grid_params1.max_y = INIParser::getInt(gridsSection, prefix + "max_y");
            grid_params1.voltage = (Const::eV / (Const::K_b * Const::EV_to_K)) * INIParser::getDouble(gridsSection, prefix + "grid_voltage");

            grids.push_back(Create_Grid(GridType::ReactConductor, domain, grid_params, grid_params1));
        }

    }
    /*std::vector<Grid*> grids;

    for (const auto& [section_name, params] : iniData) 
    {
        if (section_name.rfind("grid_", 0) == 0)  // section starts with grid_
        {
            std::string type = INIParser::getString(params, "type");

            if (type == "circular")
            {
                GridParams grid_params;
                grid_params.grid_centerx     = INIParser::getInt(params, "x_center");
                grid_params.grid_centery     = INIParser::getInt(params, "y_center");
                grid_params.num_electrodes   = INIParser::getInt(params, "electrode_number");
                grid_params.grid_rad         = INIParser::getDouble(params, "grid_radius");
                grid_params.electrode_rad    = INIParser::getDouble(params, "electrod_radius");
                grid_params.voltage          = (Const::eV/(Const::K_b*Const::EV_to_K))*INIParser::getDouble(params, "electrode_voltage");

                grids.push_back(Create_Grid(GridType::Circular, domain, grid_params));
            }
        }
    }
    */
    //---grid initialization------------------


    //-----emitter--
    std::vector<Emitter> emitters;
    const auto& emitterSection = iniData.at("Emitters");

    // Get the number of emitters
    int emitter_count = INIParser::getInt(emitterSection, "count");

    for (int i = 0; i < emitter_count; i++)
    {
        std::string prefix = "emitter_" + std::to_string(i) + ".";

        EmitterParams params;
        params.x0 = INIParser::getDouble(emitterSection, prefix + "x0");
        params.y0 = INIParser::getDouble(emitterSection, prefix + "y0");
        params.x1 = INIParser::getDouble(emitterSection, prefix + "x1");
        params.y1 = INIParser::getDouble(emitterSection, prefix + "y1");
        params.temp = INIParser::getDouble(emitterSection, prefix + "temp");
        params.numparticle = INIParser::getInt(emitterSection, prefix + "numparticle");
        params.vdx = INIParser::getDouble(emitterSection, prefix + "vdx");
        params.vdy = INIParser::getDouble(emitterSection, prefix + "vdy");
        params.species_idx1 = INIParser::getInt(emitterSection, prefix + "species_idx1");
        params.species_idx2 = INIParser::getInt(emitterSection, prefix + "species_idx2");

        emitters.emplace_back(params,domain);
    }
    ///---emitter----
    

    //initializing the species by creating instances if Init class
    for(Species &sp : species_list)
    {
        Init init(sp,domain,grids);
    } 
    
    FieldSolve fieldsolver(domain,grids);
    
    for (Species &sp:species_list)
	{
		sp.ScatterSpecies_serial();
	}

    domain.ComputeRho(species_list);
    //display::print("Rho: ", domain.rho(1,0));

    //fieldsolver.nrpcgsolver();
    //fieldsolver.pcgsolver();
    fieldsolver.PotentialSolver();
    fieldsolver.CalculateEfield();

    for (Species &sp:species_list)
	{
        sp.Rewind_species();
    }

    //species_list[0].den.display();
    //domain.rho.display();

    //int index = 0;

    //pre-calculate electron cross-section for energy level(DE_CS*1,DE_CS*2,DE_CS*3........DE_CS*(CS_RANGES-1))
    CollisionHandler testcoll(domain);
    testcoll.set_electron_cross_sections();
    //Calculate total cross-section for energy levels(DE_CS*1,DE_CS*2,DE_CS*3........DE_CS*(CS_RANGES-1))
    testcoll.calc_total_cross_sections();
    //--------------MAIN LOOP----------------------- 
    for(int ts = 0 ; ts < NUM_TS + 1; ts++)
    {
        if(ts%1==0)
        {
            for(auto &emit : emitters)
            {
                emit.inject(species_list);
            }
        }

        for (Species &sp:species_list)
		{
			sp.ScatterSpecies_serial();
		}

        domain.ComputeRho(species_list);
        //species_list[0].den.display();
        
        fieldsolver.PotentialSolver();
        fieldsolver.CalculateEfield();

        //---------particle-mover-------------------------------------------- 
        for (Species &sp:species_list)
		{
			if (sp.name == "ion" && domain.ionfixed == 1) continue;
            sp.Push_species_serial(1,grids);
		}

        //testcoll.handle_collisions(species_list[0],species_list[1]);
        //testcoll.handle_collisions(species_list[2],species_list[1]);

        if(ts%write_diagnostics == 0)
        {
            output.diagnostics(ts,species_list);
        }
        
    
        if(ts%write_interval_phase == 0)
        {
            if(ts%write_interval == 0)
            {
                
               // output.storem_to_matrix(ts,species_list);
                for(auto &sp:species_list)
                {
                    output.write_particle_data(ts,sp);
                }
            }        
        }
        if(ts%write_interval == 0)
        {
            output.storeKE_to_matrix(ts,species_list);
            if(write_flag == 1 || write_flag == 2)
            {
                output.write_field_data(ts);
                output.write_den_data(ts,species_list[0]);
                output.write_den_data(ts,species_list[1]);
            }
            
        }
    }

    output.write_ke();
    //output.write_m();
    
    auto end_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    std::cout << "Elapsed time: " << elapsed_time.count() << "seconds." <<"or "<< elapsed_time.count()/60<<"minutes"<<std::endl;

    return 0;
}
