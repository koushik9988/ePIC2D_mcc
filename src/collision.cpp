/*
Implementation file for MCC collisional models
*/

#include "collision.h"

CollisionHandler::CollisionHandler(Domain &domain, int gas, const char* data_path) : domain(domain), gas(gas)
{
    sigma_tot_e.resize(CS_RANGES, 0.0);
    sigma_tot_pion.resize(CS_RANGES, 0.0);
    sigma_tot_nion.resize(CS_RANGES, 0.0);

    sigma_ela.resize(CS_RANGES, 0.0);
    sigma_exc.resize(CS_RANGES, 0.0);
    sigma_ionz.resize(CS_RANGES, 0.0);
    sigma_mex_pion.resize(CS_RANGES, 0.0);
    sigma_mex_nion.resize(CS_RANGES, 0.0);
    sigma_e_detach.resize(CS_RANGES, 0.0);

    if (gas == HYDROGEN)
    {
        setup_hydrogen_cross_sections(data_path);
    }
    // Argon doesn't need setup (uses analytical expressions)
}

CollisionHandler::CollisionHandler(Domain &domain, int gas) : domain(domain), gas(gas)
{
    sigma_tot_e.resize(CS_RANGES, 0.0);
    sigma_tot_pion.resize(CS_RANGES, 0.0);
    sigma_tot_nion.resize(CS_RANGES, 0.0);

    sigma_ela.resize(CS_RANGES, 0.0);
    sigma_exc.resize(CS_RANGES, 0.0);
    sigma_ionz.resize(CS_RANGES, 0.0);
    sigma_mex_pion.resize(CS_RANGES, 0.0);
    sigma_mex_nion.resize(CS_RANGES, 0.0);
    sigma_e_detach.resize(CS_RANGES, 0.0);
}


void CollisionHandler::set_electron_cross_sections()
{
    std::vector<double> energy(CS_RANGES);
    energy[0] = DE_CS;
    std::generate(energy.begin() + 1, energy.end(), [i = 1]() mutable { return DE_CS * (i++); });

    if (gas == ARGON)
    {
        std::transform(energy.begin(), energy.end(), sigma_ela.begin(), [this](double energy_val) { return compute_elastic_CS_ar(energy_val, domain); });
        std::transform(energy.begin(), energy.end(), sigma_exc.begin(), [this](double energy_val) { return compute_excitation_CS_ar(energy_val, domain); });
        std::transform(energy.begin(), energy.end(), sigma_ionz.begin(), [this](double energy_val) { return compute_ionization_CS_ar(energy_val, domain); });
    }
    else if (gas == HYDROGEN)
    {
        std::transform(energy.begin(), energy.end(), sigma_ela.begin(), [this](double energy_val) { return compute_elastic_CS_h(energy_val, domain); });
        std::transform(energy.begin(), energy.end(), sigma_exc.begin(), [this](double energy_val) { return compute_excitation_CS_h(energy_val, domain); });
        std::transform(energy.begin(), energy.end(), sigma_ionz.begin(), [this](double energy_val) { return compute_ionization_CS_h(energy_val, domain); });
    }
    else
    {
        std::cerr << "Error: Unknown gas type\n";
    }
}

void CollisionHandler::set_pion_cross_sections()
{
    std::vector<double> energy(CS_RANGES);
    energy[0] = DE_CS;
    std::generate(energy.begin() + 1, energy.end(), [i = 1]() mutable { return DE_CS * (i++); });

    if (gas == HYDROGEN)
    {
        std::transform(energy.begin(), energy.end(), sigma_mex_pion.begin(), [this](double energy_val) { return compute_mex_cs_pion(energy_val, domain); });
        std::transform(energy.begin(), energy.end(), sigma_mex_nion.begin(), [this](double energy_val) { return compute_mex_cs_nion(energy_val, domain); });
        std::transform(energy.begin(), energy.end(), sigma_e_detach.begin(), [this](double energy_val) { return compute_e_detach(energy_val, domain); });
    }
}

void CollisionHandler::calc_total_cross_sections()
{
    for (size_t i = 0; i < CS_RANGES; ++i)
    {
        sigma_tot_e[i] = (sigma_ela[i] + sigma_exc[i] + sigma_ionz[i]) * domain.GAS_DENSITY;
        sigma_tot_pion[i] = sigma_mex_pion[i] * domain.GAS_DENSITY;
        sigma_tot_nion[i] = (sigma_e_detach[i] + sigma_mex_nion[i]) * domain.GAS_DENSITY;
    }
}

double CollisionHandler::max_electron_coll_freq()
{
    double e, v, nu, nu_max = 0.0;
    for (int i = 0; i < CS_RANGES; ++i)
    {
        e = i * DE_CS;
        v = sqrt(2.0 * e * Const::eV / Const::ME);
        nu = v * sigma_tot_e[i];
        if (nu > nu_max)
        {
            nu_max = nu;
        }
    }
    return nu_max;
}

void CollisionHandler::collision_electron(double xe, double ye, double &vxe, double &vye, double &vze, int eindex, Species &species1, Species &species2)
{
    const double F1 = Const::ME / (Const::ME + species2.mass);
    const double F2 = species2.mass / (Const::ME + species2.mass);
    double t0, t1, t2;
    double g, g2, gx, gy, gz, wx, wy, wz, theta, phi;
    double chi, eta, chi2, eta2, sc, cc, se, ce, st, ct, sp, cp, energy, e_sc, e_ej;
    double g_before, g_after;

    // Calculate relative velocity before collision & velocity of the center of mass before collision
    gx = vxe;
    gy = vye;
    gz = vze;
    g  = sqrt(gx * gx + gy * gy + gz * gz);
    g_before = g;
    wx = F1 * vxe;
    wy = F1 * vye;
    wz = F1 * vze;

    // Find Euler angles
    if (gx == 0)
    {
        theta = 0.5 * Const::PI;
    }
    else
    {
        theta = atan2(sqrt(gy * gy + gz * gz), gx);
    }
    if (gy == 0)
    {
        if (gz > 0)
        {
            phi = 0.5 * Const::PI;
        }
        else
        {
            phi = -0.5 * Const::PI;
        }
    }
    else
    {
        phi = atan2(gz, gy);
    }
    st = sin(theta);
    ct = cos(theta);
    sp = sin(phi);
    cp = cos(phi);

    // Choose the type of collision based on cross-sections
    t0 = sigma_ela[eindex];
    t1 = t0 + sigma_exc[eindex];
    t2 = t1 + sigma_ionz[eindex];

    if (rnd() < (t0 / t2))
    {
        // Elastic scattering
        chi = acos(1.0 - 2.0 * rnd());
        eta = 2 * Const::PI * rnd();
    }
    else if (rnd() < (t1 / t2))
    {
        // Excitation
        energy = 0.5 * Const::ME * g * g * domain.vel_norm * domain.vel_norm;
        energy = fabs(energy - E_EXC_TH * Const::eV);
        g = sqrt(2.0 * energy / Const::ME);
        g = g / domain.vel_norm;
        chi = acos(1.0 - 2.0 * rnd());
        eta = 2 * Const::PI * rnd();
    }
    else
    {
        // Ionization
        energy = 0.5 * Const::ME * g * g * domain.vel_norm * domain.vel_norm;
        energy = fabs(energy - E_ION_TH * Const::eV);
        // Shape factor w: Argon w=10 => 2w=20, Hydrogen w=8.3 => 2w=16.6
        const double shape_w = (gas == HYDROGEN) ? HYDROGEN_SHAPE_W : ARGON_SHAPE_W;
        const double shape_2w = 2.0 * shape_w;
        e_ej = 10.0 * tan(rnd() * atan(energy / Const::eV / shape_2w)) * Const::eV;
        e_sc = fabs(energy - e_ej);
        g = sqrt(2.0 * e_sc / Const::ME);
        g2 = sqrt(2.0 * e_ej / Const::ME);
        g = g / domain.vel_norm;
        g2 = g2 / domain.vel_norm;
        chi = acos(sqrt(e_sc / energy));
        chi2 = acos(sqrt(e_ej / energy));
        eta = 2 * Const::PI * rnd();
        eta2 = eta + Const::PI;

        // Compute velocity of emitted electron
        sc = sin(chi2);
        cc = cos(chi2);
        se = sin(eta2);
        ce = cos(eta2);
        gx = g2 * (ct * cc - st * sc * ce);
        gy = g2 * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
        gz = g2 * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);

        // Add the new electron and ion
        species1.AddParticle(Particle(xe, ye, (wx + F2 * gx), (wy + F2 * gy), (wz + F2 * gz)));
        double vx_new = Init::SampleVel(species2, species2.temp) / domain.vel_norm;
        double vy_new = Init::SampleVel(species2, species2.temp) / domain.vel_norm;
        double vz_new = Init::SampleVel(species2, species2.temp) / domain.vel_norm;
        species2.AddParticle(Particle(xe, ye, vx_new, vy_new, vz_new));
    }

    // Scatter the primary electron
    sc = sin(chi);
    cc = cos(chi);
    se = sin(eta);
    ce = cos(eta);

    // Compute new relative velocity
    gx = g * (ct * cc - st * sc * ce);
    gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);

    // Post-collision velocity of the colliding electron in lab frame
    vxe = wx + F2 * gx;
    vye = wy + F2 * gy;
    vze = wz + F2 * gz;

    g_after = sqrt(vxe * vxe + vye * vye + vze * vze);
    domain.delta_g = (fabs(g_after - g_before));
}


void CollisionHandler::collision_pion(double xe, double ye, double &vx_1, double &vy_1, double &vz_1, double &vx_2, double &vy_2, double &vz_2,
    int eindex, Species &species1, Species &species2)
{
    double g, gx, gy, gz, wx, wy, wz;
    double theta, phi, chi, eta, st, ct, sp, cp, sc, cc, se, ce;

    // calculate relative velocity before collision
    gx = vx_1 - vx_2;
    gy = vy_1 - vy_2;
    gz = vz_1 - vz_2;
    g  = sqrt(gx * gx + gy * gy + gz * gz);

    // in case of ion-neutral collision we have to consider neutral motion
    wx = (vx_1 * species1.mass + vx_2 * species2.mass) / (species1.mass + species2.mass);
    wy = (vy_1 * species1.mass + vy_2 * species2.mass) / (species1.mass + species2.mass);
    wz = (vz_1 * species1.mass + vz_2 * species2.mass) / (species1.mass + species2.mass);

    // find Euler angles
    if (gx == 0)
    {
        theta = 0.5 * Const::PI;
    }
    else
    {
        theta = atan2(sqrt(gy * gy + gz * gz), gx);
    }
    if (gy == 0)
    {
        if (gz > 0)
        {
            phi = 0.5 * Const::PI;
        }
        else
        {
            phi = -0.5 * Const::PI;
        }
    }
    else
    {
        phi = atan2(gz, gy);
    }

    // determine the type of collision based on cross sections and generate scattering angle
    if (rnd() < 1)
    {
        chi = acos(1.0 - 2.0 * rnd());
        eta = 2 * Const::PI * rnd();
    }

    sc = sin(chi);
    cc = cos(chi);
    se = sin(eta);
    ce = cos(eta);
    st = sin(theta);
    ct = cos(theta);
    sp = sin(phi);
    cp = cos(phi);

    // compute new relative velocity
    gx = g * (ct * cc - st * sc * ce);
    gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);

    // post-collision velocity of the ion
    double F1 = species2.mass / (species1.mass + species2.mass);
    vx_1 = wx + F1 * gx;
    vy_1 = wy + F1 * gy;
    vz_1 = wz + F1 * gz;
}


void CollisionHandler::collision_nion(double xe, double ye, double &vx_1, double &vy_1, double &vz_1, double &vx_2, double &vy_2, double &vz_2,
    int eindex, Species &species1, Species &species2, Species &electron_species)
{
    double gx = vx_1 - vx_2;
    double gy = vy_1 - vy_2;
    double gz = vz_1 - vz_2;
    double g  = sqrt(gx * gx + gy * gy + gz * gz);

    double m1 = species1.mass;
    double m2 = species2.mass;
    double m_total = m1 + m2;

    double wx = (vx_1 * m1 + vx_2 * m2) / m_total;
    double wy = (vy_1 * m1 + vy_2 * m2) / m_total;
    double wz = (vz_1 * m1 + vz_2 * m2) / m_total;

    // Cross sections for momentum exchange + detachment
    double t0 = sigma_mex_nion[eindex];
    double t1 = t0 + sigma_e_detach[eindex];
    double r  = rnd();

    if (r < (t0 / t1))
    {
        // Momentum exchange (same as collision_ion)
        double theta, phi, chi, eta;
        if (gx == 0) theta = 0.5 * Const::PI;
        else theta = atan2(sqrt(gy * gy + gz * gz), gx);

        if (gy == 0)
            phi = (gz > 0) ? 0.5 * Const::PI : -0.5 * Const::PI;
        else
            phi = atan2(gz, gy);

        chi = acos(1.0 - 2.0 * rnd());
        eta = 2.0 * Const::PI * rnd();

        double st = sin(theta), ct = cos(theta);
        double sp = sin(phi),   cp = cos(phi);
        double sc = sin(chi),   cc = cos(chi);
        double se = sin(eta),   ce = cos(eta);

        // new relative velocity
        gx = g * (ct * cc - st * sc * ce);
        gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
        gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);

        double F1 = m2 / (m1 + m2);
        vx_1 = wx + F1 * gx;
        vy_1 = wy + F1 * gy;
        vz_1 = wz + F1 * gz;
    }
    else
    {
        // Electron detachment collision
        double g_mag = sqrt(gx * gx + gy * gy + gz * gz);
        double mu = (m1 * m2) / m_total;
        double E_kin_com = 0.5 * mu * (g_mag * g_mag * domain.vel_norm * domain.vel_norm);

        double E_th_det = 2.25 * Const::eV;
        if (E_kin_com < E_th_det) return;

        double E_electron_ej = 0.0;
        double v_e_com = sqrt(2.0 * E_electron_ej / Const::ME);

        double chi = acos(1.0 - 2.0 * rnd());
        double eta = 2.0 * Const::PI * rnd();

        double sc = sin(chi), cc = cos(chi);
        double se = sin(eta), ce = cos(eta);

        double ve_com_x = v_e_com * sc * ce;
        double ve_com_y = v_e_com * sc * se;
        double ve_com_z = v_e_com * cc;

        electron_species.AddParticle(Particle(xe, ye,
            (wx + ve_com_x / domain.vel_norm),
            (wy + ve_com_y / domain.vel_norm),
            (wz + ve_com_z / domain.vel_norm)));
    }
}


void CollisionHandler::handle_collisions(Species &projectile, Species &target_gas, Species &electron_species)
{
    for (auto it = projectile.part_list.begin(); it != projectile.part_list.end(); )
    {
        Particle &part = *it;

        double v_sqr = (part.vel[0] * part.vel[0] + part.vel[1] * part.vel[1] + part.vel[2] * part.vel[2]) * domain.vel_norm * domain.vel_norm;
        double velocity = sqrt(v_sqr);
        double energy = (0.5 * projectile.mass * v_sqr) / Const::eV;
        int energy_index = std::min(static_cast<int>(energy / DE_CS + 0.5), static_cast<int>(CS_RANGES) - 1);

        double sigma_tot = 0.0;

        if (projectile.name == "electron" || projectile.name == "ebeam")
        {
            sigma_tot = sigma_tot_e[energy_index];
        }
        else if (projectile.name == "ion" && projectile.charge > 0)
        {
            sigma_tot = sigma_tot_pion[energy_index];
        }
        else if ((projectile.name == "beam" || projectile.name == "negion") && projectile.charge < 0)
        {
            sigma_tot = sigma_tot_nion[energy_index];
        }

        double nu = sigma_tot * velocity;
        double p_coll = 1 - exp(-nu * (domain.DT / domain.W));

        if (p_coll > rnd())
        {
            // Maxwellian sample for target velocity (normalized)
            double vxg = Init::SampleVel(target_gas, target_gas.temp) / domain.vel_norm;
            double vyg = Init::SampleVel(target_gas, target_gas.temp) / domain.vel_norm;
            double vzg = Init::SampleVel(target_gas, target_gas.temp) / domain.vel_norm;

            if (projectile.name == "electron" || projectile.name == "ebeam")
            {
                collision_electron(part.pos[0], part.pos[1], part.vel[0], part.vel[1], part.vel[2], energy_index, projectile, target_gas);
            }
            else if (projectile.name == "ion" && projectile.charge > 0)
            {
                collision_pion(part.pos[0], part.pos[1], part.vel[0], part.vel[1], part.vel[2], vxg, vyg, vzg, energy_index, projectile, target_gas);
            }
            else if ((projectile.name == "beam" || projectile.name == "negion") && projectile.charge < 0)
            {
                collision_nion(part.pos[0], part.pos[1], part.vel[0], part.vel[1], part.vel[2], vxg, vyg, vzg, energy_index, projectile, target_gas, electron_species);
                it = projectile.part_list.erase(it);
                continue;
            }
        }

        ++it;
    }
}


double CollisionHandler::average_collision_frequency(Species &projectile)
{
    double total_nu = 0.0;
    int num_particles = projectile.part_list.size();

    if (num_particles == 0)
    {
        return 0.0;
    }

    for (const Particle &part : projectile.part_list)
    {
        double v_sqr = (part.vel[0] * part.vel[0] + part.vel[1] * part.vel[1] + part.vel[2] * part.vel[2]) * domain.vel_norm * domain.vel_norm;
        double velocity = sqrt(v_sqr);
        double energy = (0.5 * projectile.mass * v_sqr) / Const::eV;
        int energy_index = std::min(static_cast<int>(energy / DE_CS + 0.5), static_cast<int>(CS_RANGES) - 1);

        double sigma_tot = 0.0;
        if (projectile.name == "electron" || projectile.name == "ebeam")
        {
            sigma_tot = sigma_tot_e[energy_index];
        }
        else if (projectile.name == "ion" && projectile.charge > 0)
        {
            sigma_tot = sigma_tot_pion[energy_index];
        }
        else if ((projectile.name == "beam" || projectile.name == "negion") && projectile.charge < 0)
        {
            sigma_tot = sigma_tot_nion[energy_index];
        }

        double nu = sigma_tot * velocity;
        total_nu += nu;
    }

    return total_nu / num_particles;
}
