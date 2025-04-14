#include "grid.h"

Grid::Grid(Domain& domain) : domain(domain) {}

void Grid::AddElectrode(double x, double y, double r, double voltage)
{
    electrodes.emplace_back(x, y, r, voltage, domain);
}

CustomGrid::CustomGrid(Domain& domain) : Grid(domain) {}

bool CustomGrid::IsGrid(double x, double y)
{
    for (auto &electrode : electrodes)
    {
        if (electrode.Iselectrode(x, y))
        {
            return true; // Point is inside or on any electrode
        }
    }
    return false;
}

void CustomGrid::initialize()
{
    // Implement initialization if needed
}

CircularGrid::CircularGrid(Domain& domain, int center_x, int center_y, int num_electrodes,
                          double voltage, double electrode_rad, double grid_rad)
    : Grid(domain), center_x(center_x), center_y(center_y),
      num_electrodes(num_electrodes), voltage(voltage),
      electrode_rad(electrode_rad), grid_rad(grid_rad) {}

bool CircularGrid::IsGrid(double x, double y)
{
    for ( auto &electrode : electrodes)
    {
        if (electrode.Iselectrode(x, y))
        {
            return true;
        }
    }
    return false;
}

void CircularGrid::initialize()
{
    electrodes.clear();
    for (int k = 0; k < num_electrodes; ++k)
    {
        double theta = k * (2.0 * M_PI / num_electrodes);
        double x0 = center_x + grid_rad * cos(theta);
        double y0 = center_y + grid_rad * sin(theta);
        electrodes.emplace_back(x0, y0, electrode_rad, voltage, domain);
    }
}


ReactConductorGrid::ReactConductorGrid(int min_x, int max_x, int min_y, int max_y, double voltage, Domain& domain)
    : min_x(min_x), max_x(max_x), min_y(min_y), max_y(max_y), voltage(voltage), Grid(domain) {}

bool ReactConductorGrid::IsGrid(double x, double y)
{
    return (x >= min_x && x <= max_x && y >= min_y && y <= max_y);
}

void ReactConductorGrid::initialize()
{
    int min_i = static_cast<int>(domain.XtoL(min_x));
    int max_i = static_cast<int>(domain.XtoL(max_x));
    int min_j = static_cast<int>(domain.YtoL(min_y));
    int max_j = static_cast<int>(domain.YtoL(max_y));

    for (int i = min_i; i <= max_i; ++i)
    {
        for (int j = min_j; j <= max_j; ++j)
        {
            if (i == min_i || i == max_i || j == min_j || j == max_j)
            {
                domain.rho(i, j) = voltage;
            }
        }
    }
}


Grid* Create_Grid(GridType type, Domain& domain, const GridParams& params, const RectangularGridParams &rparam)
{
    switch (type)
    {
        case GridType::Custom:
        return new CustomGrid(domain);
        case GridType::Circular:
        return new CircularGrid(domain, params.grid_centerx, params.grid_centery,params.num_electrodes,
            params.voltage,params.electrode_rad, params.grid_rad);
        case GridType::Line:
        case GridType::ReactConductor:
        return new ReactConductorGrid(rparam.min_x, rparam.max_x, rparam.min_y, rparam.max_y, rparam.voltage, domain);
            return nullptr; // Not implemented yet
        default:
            return nullptr; // Handle unknown case
    }
}