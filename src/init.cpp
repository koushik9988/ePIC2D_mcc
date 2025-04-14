#include "init.h"

Init::Init(Species &species, Domain &domain, std::vector<Grid*> &grids):
species(species), domain(domain), grids(grids)
{
    auto initialization = INIParser::loadtypeextract(species.initialization);
    auto init_type = initialization.first;
    int n = initialization.second;

    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    
    double mean_x = (domain.x0 + domain.Lx) / 2.0;
    double mean_y = (domain.y0 + domain.Ly) / 2.0;
    double sigma_x = domain.Lx / 5;
    double sigma_y = domain.Ly / 5;
    
    std::normal_distribution<> gaussian_x(mean_x, sigma_x);
    std::normal_distribution<> gaussian_y(mean_y, sigma_y);

    for (int p = 0; p < species.numparticle; p++)
    {
        double x = 0, y = 0;
        double vx = 0, vy = 0, vz = 0;

        bool particle_on_grid = false;
    
        while (true) // We keep trying until we get a valid position
        {
            // Generate particle position based on the initialization type
            if (init_type == "random")         
            {             
                x = domain.x0 + domain.Lx * rnd();             
                y = domain.y0 + domain.Ly * rnd();         
            }          
            else if (init_type == "uniform")         
            {             
                int nx = sqrt(species.numparticle);             
                int px = p % nx;             
                int py = p / nx;             
                x = domain.x0 + px * (domain.Lx / (nx - 1));  // Avoid grid boundary by adjusting range           
                y = domain.y0 + py * (domain.Ly / (nx - 1));  // Avoid grid boundary by adjusting range         
            }          

            else if (init_type == "sin")         
            {             
                x = domain.x0 + p * (domain.Lx / (species.numparticle - 1));             
                y = domain.y0 + p * (domain.Ly / (species.numparticle - 1));             
                double kx = 2 * Const::PI * n / domain.Lx;             
                double ky = 2 * Const::PI * n / domain.Ly;             
                x += sin(kx * x);             
                y += sin(ky * y);         
            }          
            else if (init_type == "cos")          
            {             
                x = domain.x0 + p * (domain.Lx / (species.numparticle - 1));             
                y = domain.y0 + p * (domain.Ly / (species.numparticle - 1));             
                double kx = 2 * Const::PI * n / domain.Lx;             
                double ky = 2 * Const::PI * n / domain.Ly;             
                x += cos(kx * x);             
                y += cos(ky * y);         
            }          

            else if (init_type == "gaussian")         
            {             
                x = gaussian_x(gen);             
                y = gaussian_y(gen);         
            }                  

            // Check if the position (x, y) is outside all the grids
            particle_on_grid = false;  // Start with assuming the particle is not inside any grid
            for (auto grid : grids)
            {
                if (grid->IsGrid(x, y)) // If the position is inside any grid
                {
                    particle_on_grid = true;  // Mark it as inside a grid
                    break;  // Exit the loop early
                }
            }

            // If the particle is not inside any grid, we can accept it
            if (!particle_on_grid) 
            {
                // Sample the velocity
                vx = Init::SampleVel(species) + species.vsx * domain.vel_norm;
                vy = Init::SampleVel(species) + species.vsy * domain.vel_norm;
                vz = 0;  // Init::SampleVel(species);  // You might want vz to be set if needed

                vx /= domain.vel_norm;
                vy /= domain.vel_norm;
                vz /= domain.vel_norm;

                // Add the particle to the species
                species.AddParticle(Particle(x, y, vx, vy, vz));
                break;  // Exit the loop since we've successfully initialized the particle
            }
        // Otherwise, the loop will continue and try generating a new position
        }
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
