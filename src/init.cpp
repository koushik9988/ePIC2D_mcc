#include "init.h"

Init::Init(Species &species, Domain &domain) : species(species), domain(domain) {
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
            x = domain.x0 + px * (domain.Ly / (nx - 1));
            y = domain.y0 + py * (domain.Ly / (nx - 1));
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
        
        vx = Init::SampleVel(species) + species.vsx * domain.vel_norm;
        vy = Init::SampleVel(species) + species.vsy * domain.vel_norm;
        vz = 0;//Init::SampleVel(species);

        vx /= domain.vel_norm;
        vy /= domain.vel_norm;
        vz /= domain.vel_norm;

        //display::print("Particle Position: ", x, ":",y);

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


// Function to inject particles
void inject(Species &species, Domain &domain, int numParticles)
{
    // Unused variable, removed to clear warning
    // double v_th = sqrt(2 * Const::K_b * domain.tempE * Const::EV_to_K / Const::ME);

    // Particle injection from the left side
    for (int p = 0; p < numParticles; p++)
    {
        
        double x = domain.x0+1;
        double y = domain.unirand(8,24);
        //double y = (rnd())*domain.Ly;
        double vx = Init::SampleVel(species);
        double vy = Init::SampleVel(species);
        double vz = 0;
        vx /= domain.vel_norm;
        vy /= domain.vel_norm;
        vz /= domain.vel_norm;
        species.AddParticle(Particle(x, y, vx, vy,vz));
    }

    // Particle injection from the right side
    for (int p = 0; p < numParticles; p++)
    {
        double x = domain.Lx -1;
        double y = domain.unirand(8,24);
    
        //double y = (rnd())*domain.Ly;

        double vx = Init::SampleVel(species);
        double vy = Init::SampleVel(species);
        double vz = 0;
        vx /= domain.vel_norm;
        vy /= domain.vel_norm;
        vx /= domain.vel_norm;
        species.AddParticle(Particle(x, y, vx, vy,vz));
    }

    // Particle injection from the bottom side
    for (int p = 0; p < numParticles; p++)
    {   
        double x = domain.unirand(8,24);
        
        //double x = (rnd())*domain.Lx;
        double y = 0+1;

        double vx = Init::SampleVel(species);
        double vy = Init::SampleVel(species);
        double vz = 0;
        vx /= domain.vel_norm;
        vy /= domain.vel_norm;
        vz /= domain.vel_norm;
        species.AddParticle(Particle(x, y, vx, vy,vz));
    }

    // Particle injection from the top side
    for (int p = 0; p < numParticles; p++)
    {
        double x = domain.unirand(8,24);
        
        //double x = (rnd())*domain.Lx;
        double y = domain.Lx - 1;
        double vx = Init::SampleVel(species);
        double vy = Init::SampleVel(species);
        double vz = 0;
        vx /= domain.vel_norm;
        vy /= domain.vel_norm;
        vz /= domain.vel_norm;
        species.AddParticle(Particle(x, y, vx, vy,vz));
    }
}
