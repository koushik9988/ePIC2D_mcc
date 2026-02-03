#include "init.h"

//Init::Init(Species &species, Domain &domain, std::vector<Grid*> &grids):


Init::Init(Species& species, Domain& domain, std::vector<Grid*>& grids) : species(species), domain(domain), grids(grids)
{
    auto initilization1 = INIParser::loadtypeextract(species.initialization_posx);
    auto initilization2 = INIParser::loadtypeextract(species.initialization_posy);
    auto initilization3 = INIParser::loadtypeextract(species.initialization_velx);
    auto initilization4 = INIParser::loadtypeextract(species.initialization_vely);


    auto [posx_type, n1, A1] = initilization1;
    auto [posy_type, n2, A2] = initilization2;
    auto [velx_type, n3, A3] = initilization3;
    auto [vely_type, n4, A4] = initilization4;

    for (int p = 0; p < species.numparticle; ++p)
    {
        double x = 0.0, y = 0.0;
        double vx = 0.0, vy = 0.0, vz = 0.0;
        bool valid_position = false;

        // --- Position initialization ---
        while (!valid_position)
        {
            
            if (posx_type == "uniform" && posy_type == "uniform")
            {
                int nx = static_cast<int>(sqrt(species.numparticle));
                if (nx * nx != species.numparticle)
                {
                    throw std::runtime_error("squre root of species particel number must be a perfect positive integer");
                }

                // Map particle index to grid coordinates
                int px = p % nx; // x-index in grid
                int py = p / nx; // y-index in grid

                // Position particles on a uniform grid
                x = domain.x0 + px * (domain.Lx / (nx)); // x-coordinate
                y = domain.y0 + py * (domain.Ly / (nx)); // y-coordinate
            }
            else
            {
                // --- x-position ---
                if (posx_type == "random")
                {
                    x = domain.x0 + domain.Lx * rnd();
                }
                else if (posx_type == "uniform")
                {
                    if (species.numparticle <= 1)
                    {
                        x = domain.x0;
                    }
                    else
                    {
                        x = domain.x0 + p * (domain.Lx / (species.numparticle - 1));
                    }
                }
                else if (posx_type == "sin")
                {
                    double x0 = domain.x0 + p * (domain.Lx / (species.numparticle - 1));
                    double kx = 2 * Const::PI * n1 / domain.Lx;
                    x = x0 + A1 * sin(kx * x0);
                }
                else if (posx_type == "cos")
                {
                    double x0 = domain.x0 + p * (domain.Lx / (species.numparticle - 1));
                    double kx = 2 * Const::PI * n1 / domain.Lx;
                    x = x0 + A1 * cos(kx * x0);
                }
                else if (posx_type == "extend")
                {
                    int start = n1;
                    double end = A1;
                    x = start + (end - start) * rnd();
                }

                // --- y-position ---
                if (posy_type == "random")
                {
                    y = domain.y0 + domain.Ly * rnd();
                }
                else if (posy_type == "uniform")
                {
                    if (species.numparticle <= 1)
                    {
                        y = domain.y0;
                    }
                    else
                    {
                        y = domain.y0 + p * (domain.Ly / (species.numparticle - 1));
                    }
                }
                else if (posy_type == "sin")
                {
                    double y0 = domain.y0 + p * (domain.Ly / (species.numparticle - 1));
                    double ky = 2 * Const::PI * n2 / domain.Ly;
                    y = y0 + A2 * sin(ky * y0);
                }
                else if (posy_type == "cos")
                {
                    double y0 = domain.y0 + p * (domain.Ly / (species.numparticle - 1));
                    double ky = 2 * Const::PI * n2 / domain.Ly;
                    y = y0 + A2 * cos(ky * y0);
                }
                else if (posy_type == "extend")
                {
                    int start = n2;
                    double end = A2;
                    y = start + (end - start) * rnd();
                }
            }


            if (x >= domain.Lx)
            {
                x -= domain.Lx;
            }
            else if (x < 0)
            {
                x += domain.Lx;
            }

            if (y >= domain.Ly)
            {
                y -= domain.Ly;
            }
            else if (y < 0)
            {
                y += domain.Ly;
            }

            // Check if the particle lies outside all the grids
            valid_position = true;
            for (auto& grid : grids)
            {
                if (grid->IsGrid(x, y))
                {
                    valid_position = false;
                    break;
                }
            }
        }

        // --- Velocity initialization ---
        vx = Init::SampleVel(species) + species.vsx * domain.vel_norm;
        vy = Init::SampleVel(species) + species.vsy * domain.vel_norm;

        if (velx_type == "sin")
        {
            double kx = 2 * Const::PI * n3 / domain.Lx;
            vx += A3 * domain.vel_norm * sin(kx * x);
        }
        else if (velx_type == "cos")
        {
            double kx = 2 * Const::PI * n3 / domain.Lx;
            vx += A3 * domain.vel_norm * cos(kx * x);
        }

        if (vely_type == "sin")
        {
            double ky = 2 * Const::PI * n4 / domain.Ly;
            vy += A4 * domain.vel_norm * sin(ky * y);
        }
        else if (vely_type == "cos")
        {
            double ky = 2 * Const::PI * n4 / domain.Ly;
            vy += A4 * domain.vel_norm * cos(ky * y);
        }

        // Normalize velocities
        vx /= domain.vel_norm;
        vy /= domain.vel_norm;
        vz = 0.0;

        species.AddParticle(Particle(x, y, vx, vy, vz));
    }
}

double Init::SampleVel(Species &species)
{
    double v_th = sqrt(2 * Const::K_b * species.temp * Const::EV_to_K / species.mass);
    return v_th * sqrt(2) * (rnd() + rnd() + rnd() - 1.5);
}

double Init::SampleVel(Species &species, double temp)
{
    double v_th = sqrt(2 * Const::K_b * temp * Const::EV_to_K / species.mass);
    return v_th * sqrt(2) * (rnd() + rnd() + rnd() - 1.5);
}
