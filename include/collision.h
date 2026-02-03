#ifndef _COLLISIONAL_MODEL_H_
#define _COLLISIONAL_MODEL_H_

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "domain.h"
#include "species.h"
#include "init.h"
#include "slap.h"

#include "Ar_cs.h"
#include "hydrogen_cs.h"

// Define gas types
#define ARGON 0
#define HYDROGEN 1

constexpr int CS_RANGES = 2000000; // CS_RANGES is a compile-time constant
constexpr double DE_CS = 0.005;
using cross_section = std::vector<double>;

/**
 * @class CollisionHandler
 * @brief Handles collision processes between particles and background gas in the simulation.
 */
class CollisionHandler
{
    private:
    Domain &domain;

    // total cross section for collision of electron
    cross_section sigma_tot_e;
    // total cross section for collision of positive ion and positive ion beam
    cross_section sigma_tot_pion;
    // total cross section for collision of negative ion and negative ion beam
    cross_section sigma_tot_nion;

    // electron elastic cross section
    cross_section sigma_ela;
    // electron excitation cross section
    cross_section sigma_exc;
    // electron ionization cross section
    cross_section sigma_ionz;
    // positive ion elastic cross section
    cross_section sigma_mex_pion;
    // negative ion elastic cross section
    cross_section sigma_mex_nion;
    // cross section for electron detachment from a negative ion
    cross_section sigma_e_detach;

    public:
    // type of background gas
    int gas;

    /**
     * @brief Constructor for the CollisionHandler class.
     * @param domain Reference to the simulation domain.
     * @param gas Type of background gas (e.g., ARGON, HYDROGEN).
     * @param data_path Path to the directory containing cross-section data files.
     */
    CollisionHandler(Domain& domain, int gas, const char* data_path);

    /**
     * @brief Overloaded Constructor for the CollisionHandler class.
     * @param domain Reference to the simulation domain.
     * @param gas Type of background gas (e.g., ARGON, HYDROGEN).
     */
    CollisionHandler(Domain &domain, int gas);

    /**
     * @brief function to set electron cross sections based on the selected gas type.
     */
    void set_electron_cross_sections();

    /**
     * @brief function to set positive ion cross sections based on the selected gas type.
     */
    void set_pion_cross_sections();

    /**
     * @brief function to calculate total cross sections from individual processes.
     */
    void calc_total_cross_sections();

    /**
     * @brief function to handle collisions between electron/electron beam and background gas.
     */
    void collision_electron(double xe, double ye, double &vxe, double &vye, double &vze, int eindex, Species &species1, Species &species2);

    /**
     * @brief function to handle collisions between positive ion/positive ion beam and background gas.
     */
    void collision_pion(double xe, double ye, double &vx_1, double &vy_1, double &vz_1, double &vx_2, double &vy_2, double &vz_2, int eindex, Species &species1, Species &species2);

    /**
     * @brief function to handle collisions between negative ion/negative ion beam and background gas.
     */
    void collision_nion(double xe, double ye, double &vx_1, double &vy_1, double &vz_1, double &vx_2, double &vy_2, double &vz_2, int eindex, Species &species1, Species &species2, Species &electron_species);

    /**
     * @brief function to process collisions for all particles in the projectile species against the background gas.
     */
    void handle_collisions(Species &projectile, Species &target_gas, Species &electron_species);

    /**
     * @brief function to compute the maximum electron collision frequency.
     */
    double max_electron_coll_freq();

    /**
     * @brief function to compute the average collision frequency for a given projectile species.
     */
    double average_collision_frequency(Species &projectile);
};

#endif  // _COLLISIONAL_MODEL_H_
