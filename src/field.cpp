#include "field.h"
#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

bool FieldSolve::gsfluidsolver()
{
    double dx2 = domain.dx * domain.dx;
    double dy2 = domain.dy * domain.dy;

    int nx = domain.nx;
    int ny = domain.ny;

    domain.phi = 0;

    double L2 = 0;
    bool converged = false;

    double coeff = (2/dx2) + (2/dy2);

    // Boundary conditions
    for(int i = 0; i < nx; i++)
    {
        domain.rho(i,0) = 0;
        domain.rho(i,ny-1) = 10000 * (Const::eV / (Const::K_b * Const::EV_to_K));
    }
    for(int i = 0; i < ny; i++)
    {
        domain.rho(0,i) = 0;
        domain.rho(nx-1,i) = 0;
    }

    for(int iter = 0; iter < domain.max_iteration; iter++)
    {
        for(int i = 1; i < nx-1; i++)
        {
            for(int j = 1; j < ny-1; j++)
            {
                double phi_new = (1/coeff) * ((domain.rho(i,j) - exp(domain.phi(i,j))) * ((domain.L*domain.L)/(domain.LDe*domain.LDe)) + 
                                             ((domain.phi(i+1,j) + domain.phi(i-1,j))/dx2) + 
                                             ((domain.phi(i,j+1) + domain.phi(i,j-1))/dy2));
                domain.phi(i,j) = domain.phi(i,j) + 1.4 * (phi_new - domain.phi(i,j)); // SOR
            }
        }

        if(iter % 50 == 0)
        {
            double sum = 0;
            for(int i = 1; i < nx-1; i++)
            {
                for(int j = 1; j < ny-1; j++)
                {
                    double R = ((domain.phi(i+1,j) - 2*domain.phi(i,j) + domain.phi(i-1,j))/dx2) + 
                               ((domain.phi(i,j+1) - 2*domain.phi(i,j) + domain.phi(i,j-1))/dy2) + 
                               (domain.rho(i,j) - exp(domain.phi(i,j))) * ((domain.L*domain.L)/(domain.LDe*domain.LDe));
                    sum += R * R;
                }
            }
            L2 = sqrt(sum / (nx * ny));
            domain.error = L2;
            if(L2 < domain.tolerance) { converged = true; break; }
        }
    }

    if(!converged) cerr << "GS failed to converge, L2=" << L2 << endl;
    return converged;
}

bool FieldSolve::gssolver()
{
    double dx2 = domain.dx * domain.dx;
    double dy2 = domain.dy * domain.dy;

    int nx = domain.nx;
    int ny = domain.ny;

    domain.phi = 0;

    //Matrix<double> rhotemp(nx, ny);
    //rhotemp = 0;

    double L2 = 0;
    bool converged = false;

    double coeff = (2/dx2) + (2/dy2);


    for(int i = 0; i < nx; i++)
    {
        domain.phi(i,0) = 0; 
        domain.phi(i,ny-1) = 0; // Top/bottom
    }
    for(int i = 0; i < ny; i++)
    {
        domain.phi(0,i) = 0; 
        domain.phi(nx-1,i) = 0*(Const::eV/(Const::K_b*Const::EV_to_K)); // Left/right, adjust units if needed
    }

    for (auto grid : grids)
    {
        grid->initialize();
    }

    //interior point exclude boundary points so we need to set phi for boudary by setting interior point we thn just need to set rho.
    for(int iter = 0; iter < domain.max_iteration; iter++)
    {
        for(int i = 1; i < nx-1; i++)
        {
            for(int j = 1; j < ny-1; j++)
            {
                double phi_new = (1/coeff) * (domain.rho(i,j) * ((domain.L*domain.L)/(domain.LDe*domain.LDe)) + ((domain.phi(i+1,j) + domain.phi(i-1,j))/dx2) + 
                                             ((domain.phi(i,j+1) + domain.phi(i,j-1))/dy2));
                domain.phi(i,j) = domain.phi(i,j) + 1.4 * (phi_new - domain.phi(i,j)); // SOR
            }
        }

        if(iter % 25 == 0)
        {
            double sum = 0;
            for(int i = 1; i < nx-1; i++)
            {
                for(int j = 1; j < ny-1; j++)
                {
                    double R = ((domain.phi(i+1,j) - 2*domain.phi(i,j) + domain.phi(i-1,j))/dx2) + 
                               ((domain.phi(i,j+1) - 2*domain.phi(i,j) + domain.phi(i,j-1))/dy2) + 
                               (domain.rho(i,j)) * ((domain.L*domain.L)/(domain.LDe*domain.LDe));
                    sum += R * R;
                }
            }
            L2 = sqrt(sum / (nx * ny));
            domain.error = L2;
            if(L2 < domain.tolerance) { converged = true; break; }
        }
    }

    if(!converged) cerr << "GS failed to converge, L2=" << L2 << endl;
    return converged;
}

void FieldSolve::nrpcgsolver()
{
    double dx2 = domain.dx * domain.dx;
    double dy2 = domain.dy * domain.dy;

    int nx = domain.nx;
    int ny = domain.ny;
    int n = nx * ny;

    Matrix<double> A(n,n), rhom(nx,ny);

    for(int i = 0; i < n; ++i)
    {
        int row = i / nx;
        int col = i % nx;
        A(i,i) = -4;
        if(col > 0) A(i, i-1) = 1;
        if(col < nx-1) A(i, i+1) = 1;
        if(row > 0) A(i, i-nx) = 1;
        if(row < ny-1) A(i, i+nx) = 1;
    }

    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            rhom(i,j) = -domain.rho(i,j) * dx2 * ((domain.L*domain.L)/(domain.LDe*domain.LDe));
        }
    }

    // Boundary conditions
    for(int i = 0; i < nx; i++)
    {
        rhom(i,0) = 0;
        rhom(i,ny-1) = 0;
    }
    for(int i = 0; i < ny; i++)
    {
        rhom(0,i) = 0;
        rhom(nx-1,i) = 0;
    }

    vec rhovec = flatten(rhom);
    vec<double> rhovecx(n);
    vec<double> x(n);

    // Newton-Raphson method
    for(int iter = 0; iter < 10; iter++)
    {
        for(int i = 0; i < n; i++)
        {
            rhovecx(i) = exp(x(i)) * ((domain.L*domain.L)/(domain.LDe*domain.LDe));
        }
        rhovec = rhovec + rhovecx;
        Matrix J = A;
        vec F = A * x - rhovec;
        for(int i = 0; i < n; i++)
        {
            J(i,i) = J(i,i) - exp(x(i));
        }
        vec<double> x_ini(n);
        vec y = cg(J, x_ini, F, 10000, 1e-4);
        x = x - y;
    }

    domain.error = 0; // Could compute a proper residual if needed
    Matrix xarr = unflatten(x, nx, ny);
    domain.phi = 0;
    domain.phi = xarr;
}

void FieldSolve::pcgsolver()
{
    double dx2 = domain.dx * domain.dx;
    double dy2 = domain.dy * domain.dy;

    int nx = domain.nx;
    int ny = domain.ny;
    int n = nx * ny;

    Matrix<double> A(n,n), rhom(nx,ny);

    for(int i = 0; i < n; ++i)
    {
        int row = i / nx;
        int col = i % nx;
        A(i,i) = -4;
        if(col > 0) A(i, i-1) = 1;
        if(col < nx-1) A(i, i+1) = 1;
        if(row > 0) A(i, i-nx) = 1;
        if(row < ny-1) A(i, i+nx) = 1;
    }

    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            rhom(i,j) = -domain.rho(i,j) * dx2 * ((domain.L*domain.L)/(domain.LDe*domain.LDe));
        }
    }

    // Boundary conditions
    for(int i = 0; i < nx; i++)
    {
        rhom(i,0) = 0;
        rhom(i,ny-1) = 0;
    }
    for(int i = 0; i < ny; i++)
    {
        rhom(0,i) = 0;
        rhom(nx-1,i) = 0;
    }

    vec rhovec = flatten(rhom);
    vec<double> x(n);
    x = pcg(A, x, rhovec, domain.max_iteration, domain.tolerance);

    vec<double> err = A * x - rhovec;
    domain.error = err.norm() / n;

    Matrix xarr = unflatten(x, nx, ny);
    domain.phi = 0;
    domain.phi = xarr;
}

void FieldSolve::PotentialSolver()
{
    if(domain.SolverType == "gs") gssolver();
    if(domain.SolverType == "iter" || domain.SolverType == "iterative" || domain.SolverType == "cg") pcgsolver();
    if(domain.SolverType == "hybrid") nrpcgsolver();
    if(domain.SolverType == "gshybrid") gsfluidsolver();
}



void FieldSolve::CalculateEfield()
{
    // Electric field along x-direction (E_x = -d(phi)/dx)
    for(int i = 1; i < domain.nx - 1; i++)
    {
        for(int j = 0; j < domain.ny; j++)
        {
            domain.efx(i,j) = -(domain.phi(i+1,j) - domain.phi(i-1,j)) / (2 * domain.dx);
        }
    }

    if(domain.bc == "open" || domain.bc == "reflective")
    {
        for(int j = 0; j < domain.ny; j++)
        {
            domain.efx(0,j) = -(domain.phi(1,j) - domain.phi(0,j)) / domain.dx;             
            domain.efx(domain.nx-1,j) = -(domain.phi(domain.nx-1,j) - domain.phi(domain.nx-2,j)) / domain.dx; 
        }
    }
    if(domain.bc == "pbc")
    {
        for(int j = 0; j < domain.ny; j++)
        {
            // Since phi(0,j) = phi(nx-1,j), use points just inside the domain for consistency
            domain.efx(0,j) = -(domain.phi(1,j) - domain.phi(domain.nx-2,j)) / (2 * domain.dx);          
            domain.efx(domain.nx-1,j) = domain.efx(0,j); 
        }
    }
    
    // Electric field along y-direction (E_y = -d(phi)/dy)
    for(int i = 0; i < domain.nx; i++)
    {
        for(int j = 1; j < domain.ny - 1; j++)
        {
            domain.efy(i,j) = -(domain.phi(i,j+1) - domain.phi(i,j-1)) / (2 * domain.dy);
        }
    }
    
    if(domain.bc == "open"  || domain.bc == "reflective")
    {
        for(int i = 0; i < domain.nx; i++)
        {
            domain.efy(i,0) = -(domain.phi(i,1) - domain.phi(i,0)) / domain.dy;                
            domain.efy(i,domain.ny-1) = -(domain.phi(i,domain.ny-1) - domain.phi(i,domain.ny-2)) / domain.dy; 
        }
    }
    if(domain.bc == "pbc")
    {
        for(int i = 0; i < domain.nx; i++)
        {
            // Since phi(i,0) = phi(i,ny-1), use points just inside the domain
            domain.efy(i,0) = -(domain.phi(i,1) - domain.phi(i,domain.ny-2)) / (2 * domain.dy);         
            domain.efy(i,domain.ny-1) = domain.efy(i,0);
        }
    }
}