#include "species.h"


/* note :
"-> " symbol means dereferencing a pointer , here we created pointers to instances of class named domain.
so to dereferenced the instances we have used the symbol "->"*/
Species::Species(string name, double mass, double charge, double spwt, double temp, 
    int numparticle, double vsx, double vsy, double fract_den , std:: string initialization, Domain &domain):domain(domain)
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
    this->initialization = initialization; 
    
    den = Matrix<double>(domain.nx,domain.ny);
}

void Species::AddParticle(Particle part)
{
    part_list.push_back(part);
}

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

        double wl = domain.LDe * domain.LDe * domain.wpe * domain.wpe;

        part.vel[0] += qm * ((domain.density * Const::QE * domain.L) / (Const::EPS_0 * domain.W * domain.vel_norm)) * part_efx * (domain.DT * sub_cycle);
        part.vel[1] += qm * ((domain.density * Const::QE * domain.L) / (Const::EPS_0 * domain.W * domain.vel_norm)) * part_efy * (domain.DT * sub_cycle);

        double x_prev1 = part.pos[0];
        double y_prev1 = part.pos[1];

        part.pos[0] += ((domain.vel_norm) / (domain.L * domain.W)) * part.vel[0] * (domain.DT * sub_cycle);
        part.pos[1] += ((domain.vel_norm) / (domain.L * domain.W)) * part.vel[1] * (domain.DT * sub_cycle);

        double x_after1 = part.pos[0];
        double y_after1 = part.pos[1];

        
        if(fabs(x_after1-x_prev1)>= domain.dx || fabs(y_after1-y_prev1)>= domain.dy)
        {
            //cout<<"particle moved more than one domain length"<<endl;
            //display::print(fabs(x_after1-x_prev1));
        }

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
    /*Multiply 0.5*mass for all particles*/
    ke(0) *= 0.5*(spwt*mass);
    ke(1) *= 0.5*(spwt*mass);
    ke(2) *= 0.5*(spwt*mass);
    
    // Calculate the total thermal energy of all the cold electrons
    double Th = (species.temp*Const::eV)*(species.spwt)*species.numparticle;

    // Normalize the kinetic energy by the total thermal energy of cold electrons    
    ke(0) = ke(0)/Th;
    ke(1) = ke(1)/Th;
    ke(2) = ke(2)/Th;
    return ke;
}

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
    /*Multiply 0.5*mass for all particles*/
    p(0) *= (spwt*mass);
    p(1) *= (spwt*mass);
    p(2) *= (spwt*mass);
    
    // Calculate the total thermal energy of all the cold electrons
    //double Thp = sqrt(2*(species.temp*Const::eV)*species.spwt*species.numparticle*species.spwt*species.numparticle);
    double Thp = sqrt(species.temp*Const::eV/Const::ME)*species.spwt*species.numparticle*Const::ME;

    // Normalize momentum by thermal velocity of normalizing species    
    p(0) = p(0)/Thp;
    p(1) = p(1)/Thp;
    p(2) = p(2)/Thp;
    return p;
}


