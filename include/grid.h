#ifndef GRID_H
#define GRID_H

#include "domain.h"
#include "electrode.h"
#include <vector>

class Electrode; // Forward declaration of Electrode class
class Domain;

enum class GridType { Custom, Circular, Line, ReactConductor};

class Grid
{
    public:
    Grid(Domain &domain);
    virtual ~Grid() = default;
    void AddElectrode(double x, double y, double r, double voltage);
    virtual bool IsGrid(double x, double y) = 0;
    virtual void initialize() = 0;

protected:
    std::vector<Electrode> electrodes;
    Domain& domain;
};

class CustomGrid : public Grid
{
    public:
    CustomGrid(Domain& domain);
    bool IsGrid(double x, double y) override;
    void initialize() override;
};

class CircularGrid : public Grid
{
    private:
    int center_x, center_y;
    int num_electrodes;
    double voltage;
    double electrode_rad;
    double grid_rad;

    public:
    CircularGrid(Domain& domain, int center_x, int center_y, int num_electrodes,
                 double voltage, double electrode_rad, double grid_rad);
    bool IsGrid(double x, double y) override;
    void initialize() override;
};



class ReactConductorGrid : public Grid
{
    private:
    int min_x, max_x, min_y, max_y;
    double voltage;

    public:
    ReactConductorGrid(int min_x, int max_x, int min_y, int max_y, double voltage, Domain& domain);
    bool IsGrid(double x, double y) override;
    void initialize() override;
};



struct GridParams
{
    int grid_centerx;
    int grid_centery;
    int num_electrodes; 
    double grid_rad;
    double electrode_rad;
    double voltage; 
};

struct RectangularGridParams
{
    int min_x;
    int max_x;
    int min_y;
    int max_y;
    double voltage;
};

Grid* Create_Grid(GridType type, Domain& domain, const GridParams& params, const RectangularGridParams &rparam);


#endif // GRID_H