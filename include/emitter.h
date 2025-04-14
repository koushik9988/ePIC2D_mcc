#include "species.h"
#include "domain.h"
#include "init.h"

class Domain;
class Species;

struct EmitterParams
{
    double x0;
    double y0;
    double x1;
    double y1;
    double temp;
    int numparticle;
    double vdx;
    double vdy;
    int species_idx1;
    int species_idx2;
};

class Emitter
{
    public:
    Emitter(const EmitterParams &params, Domain &domain);
    void inject(std::vector<Species> &species_list); 

    private:
    EmitterParams params;
    //std::vector<Species> species_list;
    Domain &domain;
};
