#include "electrode.h"
#include <cmath>

Electrode::Electrode(double x0, double y0, double r, double voltage, Domain& domain) : x0(x0), y0(y0), r(r), voltage(voltage), domain(domain)
{
    const int N = 500000; // Number of points to approximate the circle
    for (int n = 0; n < N; ++n)
    {
        double theta = 2.0 * n * M_PI / N;
        int x = static_cast<int>(x0 + r * cos(theta) + 0.5);
        int y = static_cast<int>(y0 + r * sin(theta) + 0.5);
        if (x >= 0 && x < domain.nx && y >= 0 && y < domain.ny)
        {
            domain.phi(x, y) = voltage; // Set potential at boundary of electrode
            domain.dirichlet_nodes(x,y) = true; //set dirichlet boundary
        }
    }

    // Set potential inside the electrode region to zero
    for (int i = 0; i < domain.nx; ++i)
    {
        for (int j = 0; j < domain.ny; ++j)
        {
            // Check if the point (i, j) is inside the circle
            double distance = sqrt(pow(i - x0, 2) + pow(j - y0, 2));
            if (distance < r)
            {
                domain.phi(i, j) = 0; // Set inside points to zero
                domain.dirichlet_nodes(i,j) = true; //set dirichlet boundary
            }
        }
    }
}


bool Electrode::Iselectrode(double x, double y)
{
    double dx = x - x0;
    double dy = y - y0;
    double distance_squared = dx * dx + dy * dy;
    double radius_squared = r * r;
    // True if inside or on the boundary
    return distance_squared <= radius_squared;
}