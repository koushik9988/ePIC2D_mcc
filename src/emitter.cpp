#include "emitter.h"

Emitter::Emitter(const EmitterParams &params, Domain &domain) : params(params), domain(domain){}

void Emitter::inject(std::vector<Species> &species_list)
{
    for (int p = 0; p < params.numparticle;  p++)
    {
        double x = domain.unirand(params.x0, params.x1);
        double y = domain.unirand(params.y0, params.y1);

        // First particle
        double vx1 = Init::SampleVel(species_list[params.species_idx1], params.temp) + params.vdx* domain.vel_norm;
        double vy1 = Init::SampleVel(species_list[params.species_idx1], params.temp) + params.vdy* domain.vel_norm;
        double vz1 = 0.0;

        vx1 /= domain.vel_norm;
        vy1 /= domain.vel_norm;
        vz1 /= domain.vel_norm;

        species_list[params.species_idx1].AddParticle(Particle(x, y, vx1, vy1, vz1));

        // If the two species indices are different, inject the second particle
        if (params.species_idx1 != params.species_idx2)
        {
            double vx2 = Init::SampleVel(species_list[params.species_idx2], params.temp) + params.vdx * domain.vel_norm;
            double vy2 = Init::SampleVel(species_list[params.species_idx2], params.temp) + params.vdy * domain.vel_norm;
            double vz2 = 0.0;

            vx2 /= domain.vel_norm;
            vy2 /= domain.vel_norm;
            vz2 /= domain.vel_norm;

            species_list[params.species_idx2].AddParticle(Particle(x, y, vx2, vy2, vz2));
        }
    }
}