#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <fstream>
#include <vector>
#include <filesystem>
#include <string>
#include <sstream>
#include "domain.h"
#include "species.h"
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <map>
#include <slap.h>
#include "H5Cpp.h"
#include "matplotlibcpp.h"


using namespace H5;
namespace plt = matplotlibcpp;

using namespace H5;

class Domain;
class Species;


struct PlotFlags 
{
    int phase_space;         // 1 = plot phase space (x vs vx), 0 = off
    int config_space;        // 1 = plot config space (x vs y), 0 = off
    int electric_field;      // 1 = plot electric field quiver, 0 = off
    int ke_components;       // 1 = plot kinetic energy components, 0 = off
    int total_energy;        // 1 = plot total energy plot, 0 = off
    int species_index;       // species index for phase/config plots
    int potential_field;
    int density_contour;
};


class Output 
{
    public:

    std::map<std::string, Group> particle_groups;
    std::map<std::string, Group> den_subgroups;

    Output(const std::filesystem::path& outputfolder, Domain& domain);
    void write_particle_data(int ts, Species& species);
    //void write_particle_data(H5::Group& group, int ts, Species& species);
    void write_den_data(int ts,  Species& species);
    void write_field_data(int ts);

    void diagnostics(int ts,std::vector<Species> &species_list, const PlotFlags &flags);
    int precision;
    

    void write_metadata(int NC, int NUM_TS, int write_int, int write_int_phase, double DT,double density, int save_fig, int normscheme, 
        int subcycleint,double LDe, double LDi, double wpe, double wpi,int spno, double GAS_DENSITY);
    void write_species_metadata(std::vector<Species>& species_list);


    void write_ke();
    void write_m();
    void storeKE_to_matrix(int ts, std::vector<Species> &species_list);
    void storem_to_matrix(int ts, std::vector<Species> &species_list);

    Matrix<double> store_ke;
    Matrix<double> store_m;

    int sp_no ;//= species_list.size();
    int t_step;// = int(domain.NUM_TS/domain.write_interval) + 1 ;



    std::vector<double> time_steps;
    std::vector<double> kinetic_energy;
    std::vector<double> potential_energy;
    std::vector<double> total_energy;
    std::vector<double>Ke_x;
    std::vector<double>Ke_y;
    std::vector<double>Ke_z;

    private:
    std::filesystem::path outputfolder;
    Domain& domain;
    H5File file; // Declare H5::H5File object to handle HDF5 file operations
    Group field_data_group;
    Group KE_group;
    Group metadata_group;
    Group time_group;
    Group metadata_species;
   
};

namespace display
{
    template<typename T>
    void print(const T& value) 
    {
        std::cout << value << std::endl;
    }
    // Recursive template function to print multiple arguments
    template<typename T, typename... Args>
    void print(const T& value, Args&&... args) 
    {
        std::cout << value;
        print(std::forward<Args>(args)...); // Recursive call with the remaining arguments
    }
}

#endif
