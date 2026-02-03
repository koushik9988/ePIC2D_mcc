#include "species.h"


/* note :
"-> " symbol means dereferencing a pointer , here we created pointers to instances of class named domain.
so to dereferenced the instances we have used the symbol "->"*/
Species::Species(string name, double mass, double charge, double spwt, double temp, 
    int numparticle, double vsx, double vsy, double fract_den , std:: string initialization_posx, std:: string initialization_posy, 
    std:: string initialization_velx, std:: string initialization_vely, Domain &domain):domain(domain)
{
    this->name = name;
    this->mass = mass;
    this->charge = charge;
    this->spwt = spwt;
    this->temp = temp;
    this->numparticle = numparticle;
    this->vsx = vsx;
    this->vsy = vsy;
    this->fract_den = fract_den;
    this->initialization_posx = initialization_posx;
    this->initialization_posy = initialization_posy;
    this->initialization_velx = initialization_velx;
    this->initialization_vely = initialization_vely; 
    
    den = Matrix<double>(domain.nx,domain.ny);
}

void Species::AddParticle(Particle part)
{
    part_list.push_back(part);
}

/*
void Species::Push_species_serial(int sub_cycle, const std::vector<Grid*> &grids)
{
    for (auto it = part_list.begin(); it != part_list.end();)
    {
        Particle &part = *it;
        double qm = charge / mass;
       

        double lx = domain.XtoL(part.pos[0]);
        double ly = domain.YtoL(part.pos[1]);

        double part_efx = domain.Interpolate(lx, ly, domain.efx);
        double part_efy = domain.Interpolate(lx, ly, domain.efy);

        double dt = domain.DT * sub_cycle;
        
        // --- Half-step acceleration by E field
        double vx_minus = part.vel[0] + 0.5 * part_efx * dt * (qm * ((domain.density * Const::QE * domain.L) / (Const::EPS_0 * domain.W * domain.vel_norm)));
        double vy_minus = part.vel[1] + 0.5 * part_efy * dt * (qm * ((domain.density * Const::QE * domain.L) / (Const::EPS_0 * domain.W * domain.vel_norm)));
        double vz_minus = part.vel[2];

        double vx_plus, vy_plus, vz_plus;

        if (domain.B != 0.0 )
        {
            // Unnormalize velocities
            vx_minus *= domain.vel_norm;
            vy_minus *= domain.vel_norm;
            vz_minus *= domain.vel_norm;

            double B = domain.B;
            double theta = domain.theta;
            double azimuth = domain.azimuth;

            double Bx = B * sin(theta) * cos(azimuth);
            double By = B * sin(theta) * sin(azimuth); 
            double Bz = B * cos(theta);


            // t = (q/m) * B * dt / 2
            double tx = 0.5 * qm * Bx * (dt / domain.W);
            double ty = 0.5 * qm * By * (dt / domain.W);
            double tz = 0.5 * qm * Bz * (dt / domain.W);

            // v_prime = v_minus + v_minus x t
            double v_prime_x = vx_minus + (vy_minus * tz - vz_minus * ty);
            double v_prime_y = vy_minus + (vz_minus * tx - vx_minus * tz);
            double v_prime_z = vz_minus + (vx_minus * ty - vy_minus * tx);

            // s = 2 / (1 + |t|^2)
            double t2 = tx * tx + ty * ty + tz * tz;
            double s = 2.0 / (1.0 + t2);

            // v_plus = v_minus + v_prime x (s * t)
            // sx = s * tx  | sy = s * ty |  sz = s * tz
            vx_plus = vx_minus + (v_prime_y * (s * tz) - v_prime_z * (s * ty));
            vy_plus = vy_minus + (v_prime_z * (s * tx) - v_prime_x * (s * tz));
            vz_plus = vz_minus + (v_prime_x * (s * ty) - v_prime_y * (s * tx));

            // Normalize velocities
            vx_plus /= domain.vel_norm;
            vy_plus /= domain.vel_norm;
            vz_plus /= domain.vel_norm;
        }
        else
        {
            vx_plus = vx_minus;
            vy_plus = vy_minus;
            vz_plus = vz_minus;
        }

        // --- Half-step acceleration by E field (again)
        part.vel[0] = vx_plus + 0.5 * part_efx * dt * (qm * ((domain.density * Const::QE * domain.L) / (Const::EPS_0 * domain.W * domain.vel_norm)));
        part.vel[1] = vy_plus + 0.5 * part_efy * dt * (qm * ((domain.density * Const::QE * domain.L) / (Const::EPS_0 * domain.W * domain.vel_norm)));
        part.vel[2] = vz_plus;

        // --- Position update
        part.pos[0] += ((domain.vel_norm) / (domain.L * domain.W)) * part.vel[0] * dt;
        part.pos[1] += ((domain.vel_norm) / (domain.L * domain.W)) * part.vel[1] * dt;
        
        if (domain.bc == "pbc")
        {
            if (part.pos[0] < domain.x0)
            {
                part.pos[0] = part.pos[0] + domain.Lx;
            }
            else if (part.pos[0] >= domain.x0 + domain.Lx)
            {
                part.pos[0] = part.pos[0] - domain.Lx;
            }

            if (part.pos[1] < domain.x0)
            {
                part.pos[1] = part.pos[1] + domain.Ly;
            }
            else if (part.pos[1] >= domain.x0 + domain.Ly)
            {
                part.pos[1] = part.pos[1] - domain.Ly;
            }

            ////
            bool erase_particle = false;
            for (auto grid : grids)
            {
                if(grid->IsGrid(part.pos[0],part.pos[1]))
                {
                    erase_particle = true;
                    break;
                }
            }
            if (erase_particle)
            {
                it = part_list.erase(it);
                continue;
            } 
            ////
        }
        else if (domain.bc == "open")
        {
            // Remove particles that move outside the domain
            if (part.pos[0] < domain.x0 || part.pos[0] >= domain.x0 + domain.Lx ||
                part.pos[1] < domain.x0 || part.pos[1] >= domain.x0 + domain.Ly)
            {
                it = part_list.erase(it); // Remove particle and update iterator
                continue; // Skip incrementing the iterator
            }

            bool erase_particle = false;
            for (auto grid : grids)
            {
                if(grid->IsGrid(part.pos[0],part.pos[1]))
                {
                    erase_particle = true;
                    break;
                }
            }
            if (erase_particle)
            {
                it = part_list.erase(it);
                continue;
            }     
        }
        else if (domain.bc == "reflective")
        {
            // Reflect particles at the boundaries
            if (part.pos[0] < domain.x0)
            {
                part.pos[0] = domain.x0 + (domain.x0 - part.pos[0]);
                part.vel[0] = -part.vel[0]; // Reflect velocity
            }
            else if (part.pos[0] >= domain.x0 + domain.Lx)
            {
                part.pos[0] = domain.x0 + domain.Lx - (part.pos[0] - (domain.x0 + domain.Lx));
                part.vel[0] = -part.vel[0]; // Reflect velocity
            }

            if (part.pos[1] < domain.x0)
            {
                part.pos[1] = domain.x0 + (domain.x0 - part.pos[1]);
                part.vel[1] = -part.vel[1]; // Reflect velocity
            }
            else if (part.pos[1] >= domain.x0 + domain.Ly)
            {
                part.pos[1] = domain.x0 + domain.Ly - (part.pos[1] - (domain.x0 + domain.Ly));
                part.vel[1] = -part.vel[1]; // Reflect velocity
            }

            ////
            bool erase_particle = false;
            for (auto grid : grids)
            {
                if(grid->IsGrid(part.pos[0],part.pos[1]))
                {
                    erase_particle = true;
                    break;
                }
            }
            if (erase_particle)
            {
                it = part_list.erase(it);
                continue;
            } 
            ////

        }
        it++; // Increment iterator
    }
}
*/

void Species::Push_species_serial(int sub_cycle, const std::vector<Grid*> &grids)
{
    int n = static_cast<int>(part_list.size());
    if (n == 0)
    {
        return;
    }

    std::vector<char> keep(static_cast<size_t>(n), 1);

    double qm = charge / mass;
    double dt = domain.DT * sub_cycle;
    double accel_scale = qm * ((domain.density * Const::QE * domain.L) / (Const::EPS_0 * domain.W * domain.vel_norm));
    double half_accel = 0.5 * dt * accel_scale;
    double pos_scale = (domain.vel_norm) / (domain.L * domain.W);

    bool use_B = (domain.B != 0.0);
    double tx = 0.0;
    double ty = 0.0;
    double tz = 0.0;
    double s = 0.0;
    if (use_B)
    {
        double B = domain.B;
        double theta = domain.theta;
        double azimuth = domain.azimuth;

        double Bx = B * sin(theta) * cos(azimuth);
        double By = B * sin(theta) * sin(azimuth);
        double Bz = B * cos(theta);

        // t = (q/m) * B * dt / 2
        tx = 0.5 * qm * Bx * (dt / domain.W);
        ty = 0.5 * qm * By * (dt / domain.W);
        tz = 0.5 * qm * Bz * (dt / domain.W);

        // s = 2 / (1 + |t|^2)
        double t2 = tx * tx + ty * ty + tz * tz;
        s = 2.0 / (1.0 + t2);
    }

    bool bc_pbc = (domain.bc == "pbc");
    bool bc_open = (domain.bc == "open");
    bool bc_reflective = (domain.bc == "reflective");

    #pragma omp parallel for schedule(guided, 256)
    for (int idx = 0; idx < n; idx++)
    {
        Particle &part = part_list[idx];

        double lx = domain.XtoL(part.pos[0]);
        double ly = domain.YtoL(part.pos[1]);

        double part_efx = domain.Interpolate(lx, ly, domain.efx);
        double part_efy = domain.Interpolate(lx, ly, domain.efy);

        // --- Half-step acceleration by E field
        double vx_minus = part.vel[0] + half_accel * part_efx;
        double vy_minus = part.vel[1] + half_accel * part_efy;
        double vz_minus = part.vel[2];

        double vx_plus;
        double vy_plus;
        double vz_plus;

        if (use_B)
        {
            // Unnormalize velocities
            vx_minus *= domain.vel_norm;
            vy_minus *= domain.vel_norm;
            vz_minus *= domain.vel_norm;

            // v_prime = v_minus + v_minus x t
            double v_prime_x = vx_minus + (vy_minus * tz - vz_minus * ty);
            double v_prime_y = vy_minus + (vz_minus * tx - vx_minus * tz);
            double v_prime_z = vz_minus + (vx_minus * ty - vy_minus * tx);

            // v_plus = v_minus + v_prime x (s * t)
            // sx = s * tx  | sy = s * ty |  sz = s * tz
            vx_plus = vx_minus + (v_prime_y * (s * tz) - v_prime_z * (s * ty));
            vy_plus = vy_minus + (v_prime_z * (s * tx) - v_prime_x * (s * tz));
            vz_plus = vz_minus + (v_prime_x * (s * ty) - v_prime_y * (s * tx));

            // Normalize velocities
            vx_plus /= domain.vel_norm;
            vy_plus /= domain.vel_norm;
            vz_plus /= domain.vel_norm;
        }
        else
        {
            vx_plus = vx_minus;
            vy_plus = vy_minus;
            vz_plus = vz_minus;
        }

        // --- Half-step acceleration by E field (again)
        part.vel[0] = vx_plus + half_accel * part_efx;
        part.vel[1] = vy_plus + half_accel * part_efy;
        part.vel[2] = vz_plus;

        // --- Position update
        part.pos[0] += pos_scale * part.vel[0] * dt;
        part.pos[1] += pos_scale * part.vel[1] * dt;

        bool erase_particle = false;

        if (bc_pbc)
        {
            if (part.pos[0] < domain.x0)
            {
                part.pos[0] = part.pos[0] + domain.Lx;
            }
            else if (part.pos[0] >= domain.x0 + domain.Lx)
            {
                part.pos[0] = part.pos[0] - domain.Lx;
            }

            if (part.pos[1] < domain.x0)
            {
                part.pos[1] = part.pos[1] + domain.Ly;
            }
            else if (part.pos[1] >= domain.x0 + domain.Ly)
            {
                part.pos[1] = part.pos[1] - domain.Ly;
            }

            for (auto grid : grids)
            {
                if(grid->IsGrid(part.pos[0],part.pos[1]))
                {
                    erase_particle = true;
                    break;
                }
            }
        }
        else if (bc_open)
        {
            // Remove particles that move outside the domain
            if (part.pos[0] < domain.x0 || part.pos[0] >= domain.x0 + domain.Lx ||
                part.pos[1] < domain.x0 || part.pos[1] >= domain.x0 + domain.Ly)
            {
                erase_particle = true;
            }
            else
            {
                for (auto grid : grids)
                {
                    if(grid->IsGrid(part.pos[0],part.pos[1]))
                    {
                        erase_particle = true;
                        break;
                    }
                }
            }
        }
        else if (bc_reflective)
        {
            // Reflect particles at the boundaries
            if (part.pos[0] < domain.x0)
            {
                part.pos[0] = domain.x0 + (domain.x0 - part.pos[0]);
                part.vel[0] = -part.vel[0]; // Reflect velocity
            }
            else if (part.pos[0] >= domain.x0 + domain.Lx)
            {
                part.pos[0] = domain.x0 + domain.Lx - (part.pos[0] - (domain.x0 + domain.Lx));
                part.vel[0] = -part.vel[0]; // Reflect velocity
            }

            if (part.pos[1] < domain.x0)
            {
                part.pos[1] = domain.x0 + (domain.x0 - part.pos[1]);
                part.vel[1] = -part.vel[1]; // Reflect velocity
            }
            else if (part.pos[1] >= domain.x0 + domain.Ly)
            {
                part.pos[1] = domain.x0 + domain.Ly - (part.pos[1] - (domain.x0 + domain.Ly));
                part.vel[1] = -part.vel[1]; // Reflect velocity
            }

            for (auto grid : grids)
            {
                if(grid->IsGrid(part.pos[0],part.pos[1]))
                {
                    erase_particle = true;
                    break;
                }
            }
        }

        if (erase_particle)
        {
            keep[static_cast<size_t>(idx)] = 0;
        }
    }

    int write_idx = 0;
    for (int read_idx = 0; read_idx < n; read_idx++)
    {
        if (keep[static_cast<size_t>(read_idx)])
        {
            if (write_idx != read_idx)
            {
                part_list[write_idx] = part_list[read_idx];
            }
            write_idx++;
        }
    }
    if (write_idx < n)
    {
        part_list.erase(part_list.begin() + write_idx, part_list.end());
    }
}

void Species::ScatterSpecies_serial()
{
    den = 0.0;   
    // Scatter particles and accumulate density
    for(Particle& part : part_list)
    {
        double lx = domain.XtoL(part.pos[0]);  // Convert x-position to logical coordinates
        double ly = domain.YtoL(part.pos[1]);  // Convert y-position to logical coordinates
        domain.Scatter(lx, ly, spwt, den);    // Scatter particle with specific weight
        //display::print("Particle Position: ", lx, ":", ly);
        //display::print(den(1,1));
    }

    //den.display();

    // Normalize density
    double cell_area = domain.dx * domain.dy;
    double L2 = domain.L * domain.L;
    den /= (cell_area * L2 * domain.density);

    // Apply boundary conditions
    if(domain.bc == "pbc")  // Periodic boundary conditions
    {
        for(int i = 0; i < domain.nx; i++)
        {
            den(i, 0) += den(i, domain.ny-1);          // Add opposite boundary
            den(i, domain.ny-1) = den(i, 0);          // Mirror the value
        }
        for(int j = 0; j < domain.ny; j++)
        {
            den(0, j) += den(domain.nx-1, j);
            den(domain.nx-1, j) = den(0, j);
        }
    }
    else if(domain.bc == "open" || domain.bc == "reflective")  // Non-periodic boundaries
    {
        //Handle corners (multiply by 0.25)
        den(0, 0) *= 0.25;
        den(0, domain.ny-1) *= 0.25;
        den(domain.nx-1, 0) *= 0.25;
        den(domain.nx-1, domain.ny-1) *= 0.25;

        //Handle edges (multiply by 0.5, excluding corners)
        for(int i = 1; i < domain.nx-1; i++)
        {
            den(i, 0) *= 0.5;
            den(i, domain.ny-1) *= 0.5;
        }
        for(int j = 1; j < domain.ny-1; j++)
        {
            den(0, j) *= 0.5;
            den(domain.nx-1, j) *= 0.5;
        }
    }
}


/*
void Species::Rewind_species()
{
    for (Particle &part: part_list)
    {
        double qm = charge/mass;
        double lx = domain.XtoL(part.pos[0]);
        double ly = domain.YtoL(part.pos[1]);

        double part_efx = domain.Interpolate(lx, ly, domain.efx);
        double part_efy = domain.Interpolate(lx, ly, domain.efy);

        double wl = domain.LDe*domain.LDe*domain.wpe*domain.wpe;

        part.vel[0] -= 0.5*qm*((domain.density*Const::QE*domain.L)/(Const::EPS_0*domain.W*domain.vel_norm))*part_efx*domain.DT;
        part.vel[1] -= 0.5*qm*((domain.density*Const::QE*domain.L)/(Const::EPS_0*domain.W*domain.vel_norm))*part_efy*domain.DT;
        //part.vel[2] -= 0.5*qm*((domain.density*Const::QE*domain.L)/(Const::EPS_0*domain.W*domain.vel_norm))*part_efy*domain.DT;         
    }
}
*/

void Species::Rewind_species()
{
    int n = static_cast<int>(part_list.size());
    double qm = charge/mass;
    double accel = 0.5 * qm * ((domain.density*Const::QE*domain.L)/(Const::EPS_0*domain.W*domain.vel_norm)) * domain.DT;

    #pragma omp parallel for schedule(guided, 256)
    for (int i = 0; i < n; i++)
    {
        Particle &part = part_list[i];
        double lx = domain.XtoL(part.pos[0]);
        double ly = domain.YtoL(part.pos[1]);

        double part_efx = domain.Interpolate(lx, ly, domain.efx);
        double part_efy = domain.Interpolate(lx, ly, domain.efy);

        part.vel[0] -= accel * part_efx;
        part.vel[1] -= accel * part_efy;
        //part.vel[2] -= accel * part_efy;         
    }
}

/*
vec<double> Species::Compute_KE(Species &species)
{
    vec<double> ke(3);
    for (Particle &part:part_list)
    {
        // un-normalize the velocity by multiplying with the cold thermal velocity
        ke(0) += (part.vel[0]*part.vel[0])*(domain.vel_norm)*(domain.vel_norm);
        ke(1) += (part.vel[1]*part.vel[1])*(domain.vel_norm)*(domain.vel_norm);
        ke(2) += (part.vel[2]*part.vel[2])*(domain.vel_norm)*(domain.vel_norm);
    }
    //Multiply 0.5*mass for all particles
    ke(0) *= 0.5*(spwt*mass);
    ke(1) *= 0.5*(spwt*mass);
    ke(2) *= 0.5*(spwt*mass);
    
    // Calculate the total thermal energy of all the cold electrons
    double Th = (species.temp*Const::eV)*(species.spwt)*species.numparticle;
    if(domain.normscheme == 5)
    {Th = 1;}

    // Normalize the kinetic energy by the total thermal energy of cold electrons    
    ke(0) = ke(0)/Th;
    ke(1) = ke(1)/Th;
    ke(2) = ke(2)/Th;
    return ke;
}
*/

vec<double> Species::Compute_KE(Species &species)
{
    vec<double> ke(3);
    int n = static_cast<int>(part_list.size());
    double vnorm2 = domain.vel_norm * domain.vel_norm;

    double sum0 = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;

    #pragma omp parallel for reduction(+:sum0,sum1,sum2) schedule(guided, 256)
    for (int i = 0; i < n; i++)
    {
        const Particle &part = part_list[i];
        sum0 += (part.vel[0] * part.vel[0]) * vnorm2;
        sum1 += (part.vel[1] * part.vel[1]) * vnorm2;
        sum2 += (part.vel[2] * part.vel[2]) * vnorm2;
    }

    ke(0) = sum0;
    ke(1) = sum1;
    ke(2) = sum2;

    /*Multiply 0.5*mass for all particles*/
    ke(0) *= 0.5 * (spwt * mass);
    ke(1) *= 0.5 * (spwt * mass);
    ke(2) *= 0.5 * (spwt * mass);
    
    // Calculate the total thermal energy of all the cold electrons
    double Th = (species.temp * Const::eV) * (species.spwt) * species.numparticle;
    if(domain.normscheme == 5)
    {Th = 1;}

    // Normalize the kinetic energy by the total thermal energy of cold electrons    
    ke(0) = ke(0)/Th;
    ke(1) = ke(1)/Th;
    ke(2) = ke(2)/Th;
    return ke;
}


/*
vec<double> Species::Compute_Momentum(Species &species)
{
    vec<double> p(3);
    for (Particle &part:part_list)
    {
        // un-normalize the velocity by multiplying with the cold thermal velocity
        p(0) += (part.vel[0])*domain.vel_norm;
        p(1) += (part.vel[1])*domain.vel_norm;
        p(2) += (part.vel[2])*domain.vel_norm;
    }
    //Multiply 0.5*mass for all particles
    p(0) *= (spwt*mass);
    p(1) *= (spwt*mass);
    p(2) *= (spwt*mass);
    
    // Calculate the total thermal energy of all the cold electrons
    //double Thp = sqrt(2*(species.temp*Const::eV)*species.spwt*species.numparticle*species.spwt*species.numparticle);
    double Thp = sqrt(species.temp*Const::eV/Const::ME)*species.spwt*species.numparticle*species.mass;

    if(domain.normscheme == 5)
    {Thp = 1;}
    
    // Normalize momentum by thermal velocity of normalizing species    
    p(0) = p(0)/Thp;
    p(1) = p(1)/Thp;
    p(2) = p(2)/Thp;
    return p;
}
*/

vec<double> Species::Compute_Momentum(Species &species)
{
    vec<double> p(3);
    int n = static_cast<int>(part_list.size());

    double sum0 = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;

    #pragma omp parallel for reduction(+:sum0,sum1,sum2) schedule(guided, 256)
    for (int i = 0; i < n; i++)
    {
        const Particle &part = part_list[i];
        sum0 += part.vel[0] * domain.vel_norm;
        sum1 += part.vel[1] * domain.vel_norm;
        sum2 += part.vel[2] * domain.vel_norm;
    }

    p(0) = sum0;
    p(1) = sum1;
    p(2) = sum2;

    //Multiply 0.5*mass for all particles
    p(0) *= (spwt * mass);
    p(1) *= (spwt * mass);
    p(2) *= (spwt * mass);
    
    // Calculate the total thermal energy of all the cold electrons
    //double Thp = sqrt(2*(species.temp*Const::eV)*species.spwt*species.numparticle*species.spwt*species.numparticle);
    double Thp = sqrt(species.temp*Const::eV/Const::ME)*species.spwt*species.numparticle*species.mass;

    if(domain.normscheme == 5)
    {Thp = 1;}
    
    // Normalize momentum by thermal velocity of normalizing species    
    p(0) = p(0)/Thp;
    p(1) = p(1)/Thp;
    p(2) = p(2)/Thp;
    return p;
}


