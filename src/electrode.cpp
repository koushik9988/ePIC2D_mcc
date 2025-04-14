#include "electrode.h"
#include <cmath>

Electrode::Electrode(double x0, double y0, double r, double voltage, Domain& domain) : x0(x0), y0(y0), r(r), voltage(voltage), domain(domain)
{
    const int N = 100; // Number of points to approximate the circle
    for (int n = 0; n < N; ++n)
    {
        double theta = 2.0 * n * M_PI / N;
        int x = static_cast<int>(x0 + r * cos(theta) + 0.5);
        int y = static_cast<int>(y0 + r * sin(theta) + 0.5);
        if (x >= 0 && x < domain.nx && y >= 0 && y < domain.ny)
        {
            domain.rho(x, y) = voltage; // Set potential at boundary of electrode
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