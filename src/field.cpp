#include "field.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fftw3.h>


using namespace std;

void FieldSolve::spectral()
{
    int Nx = domain.nx;
    int Ny = domain.ny;
    int Nh = (Ny / 2) + 1;
    double norm = 1.0 / (Nx * Ny);
    double Lx = domain.Lx;
    double Ly = domain.Ly;

    double dx = domain.dx;
    double dy = domain.dy;

    // Allocate real-space and Fourier-space arrays
    double* rho = fftw_alloc_real(Nx * Ny);
    double* phi = fftw_alloc_real(Nx * Ny);
    fftw_complex* rho_k = fftw_alloc_complex(Nx * Nh);
    fftw_complex* phi_k = fftw_alloc_complex(Nx * Nh);

    // Zero out all memory
    std::memset(rho, 0, sizeof(double) * Nx * Ny);
    std::memset(phi, 0, sizeof(double) * Nx * Ny);
    std::memset(rho_k, 0, sizeof(fftw_complex) * Nx * Nh);
    std::memset(phi_k, 0, sizeof(fftw_complex) * Nx * Nh);

    // Fill rho and apply normalization
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            rho[i * Ny + j] =  -domain.rho(i,j) * ((domain.L*domain.L)/(domain.LDe*domain.LDe));
        }
    }

    // Create FFTW plans
    fftw_plan forward = fftw_plan_dft_r2c_2d(Nx, Ny, rho, rho_k, FFTW_ESTIMATE);
    fftw_plan backward = fftw_plan_dft_c2r_2d(Nx, Ny, phi_k, phi, FFTW_ESTIMATE);

    // Forward transform
    fftw_execute(forward);

    // Solve Poisson in Fourier space
    for (int i = 0; i < Nx; ++i)
    {
        int kx_index = (i <= Nx / 2) ? i : i - Nx;
        double kx = 2.0 * M_PI * kx_index / Lx;

        for (int j = 0; j < Nh; ++j)
        {
            double ky = 2.0 * M_PI * j / Ly;
            //double k2 = kx * kx + ky * ky; // 
            double k2 = 4.0 / (dx * dx) * pow(sin(0.5 * kx * dx), 2) + 4.0 / (dy * dy) * pow(sin(0.5 * ky * dy), 2);

            int idx = i * Nh + j;
            if (k2 == 0)
            {
                phi_k[idx][0] = 0.0;
                phi_k[idx][1] = 0.0;
            }
            else
            {
                phi_k[idx][0] = rho_k[idx][0] / (-k2);
                phi_k[idx][1] = rho_k[idx][1] / (-k2);
            }
        }
    }

    // Inverse FFT
    fftw_execute(backward);

    domain.phi = 0;

    // Normalize and store back
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            domain.phi(i, j) = phi[i * Ny + j] * norm;
        }
    }

    // Cleanup
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    fftw_free(rho_k);
    fftw_free(phi_k);
    fftw_free(rho);
    fftw_free(phi);
    fftw_cleanup();
}



/*
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
        domain.phi(nx-1,i) = 0*domain.volatge_norm; // Left/right, adjust units if needed
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
                if(domain.dirichlet_nodes(i,j)) continue;
                double phi_new = (1/coeff) * (domain.rho(i,j) * ((domain.L*domain.L)/(domain.LDe*domain.LDe)) 
                + ((domain.phi(i+1,j) + domain.phi(i-1,j))/dx2) + ((domain.phi(i,j+1) + domain.phi(i,j-1))/dy2));
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
                    if(domain.dirichlet_nodes(i,j)) continue;
                    double R = ((domain.phi(i+1,j) - 2*domain.phi(i,j) + domain.phi(i-1,j))/dx2) + ((domain.phi(i,j+1) - 2*domain.phi(i,j) + domain.phi(i,j-1))/dy2) + 
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
*/

bool FieldSolve::gssolver()
{
    double dx2 = domain.dx * domain.dx;
    double dy2 = domain.dy * domain.dy;

    int nx = domain.nx;
    int ny = domain.ny;

    domain.phi = 0;

    double L2 = 0;
    bool converged = false;

    double coeff = (2/dx2) + (2/dy2);
    double rho_scale = (domain.L * domain.L) / (domain.LDe * domain.LDe);
    double omega = 1.4;

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < nx; i++)
    {
        domain.phi(i,0) = 0; 
        domain.phi(i,ny-1) = 0; // Top/bottom
    }
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < ny; i++)
    {
        domain.phi(0,i) = 0; 
        domain.phi(nx-1,i) = 0*domain.volatge_norm; // Left/right, adjust units if needed
    }

    for (auto grid : grids)
    {
        grid->initialize();
    }

    //interior point exclude boundary points so we need to set phi for boudary by setting interior point we thn just need to set rho.
    for(int iter = 0; iter < domain.max_iteration; iter++)
    {
        // Red sweep (i + j even)
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < nx-1; i++)
        {
            for(int j = 1; j < ny-1; j++)
            {
                if(((i + j) & 1) != 0) continue;
                if(domain.dirichlet_nodes(i,j)) continue;
                double phi_new = (1/coeff) * (domain.rho(i,j) * rho_scale 
                + ((domain.phi(i+1,j) + domain.phi(i-1,j))/dx2) + ((domain.phi(i,j+1) + domain.phi(i,j-1))/dy2));
                domain.phi(i,j) = domain.phi(i,j) + omega * (phi_new - domain.phi(i,j)); // SOR
            }
        }

        // Black sweep (i + j odd)
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < nx-1; i++)
        {
            for(int j = 1; j < ny-1; j++)
            {
                if(((i + j) & 1) == 0) continue;
                if(domain.dirichlet_nodes(i,j)) continue;
                double phi_new = (1/coeff) * (domain.rho(i,j) * rho_scale 
                + ((domain.phi(i+1,j) + domain.phi(i-1,j))/dx2) + ((domain.phi(i,j+1) + domain.phi(i,j-1))/dy2));
                domain.phi(i,j) = domain.phi(i,j) + omega * (phi_new - domain.phi(i,j)); // SOR
            }
        }

        if(iter % 25 == 0)
        {
            double sum = 0;
            #pragma omp parallel for collapse(2) reduction(+:sum) schedule(static)
            for(int i = 1; i < nx-1; i++)
            {
                for(int j = 1; j < ny-1; j++)
                {
                    if(domain.dirichlet_nodes(i,j)) continue;
                    double R = ((domain.phi(i+1,j) - 2*domain.phi(i,j) + domain.phi(i-1,j))/dx2) + ((domain.phi(i,j+1) - 2*domain.phi(i,j) + domain.phi(i,j-1))/dy2) + 
                    (domain.rho(i,j)) * rho_scale;
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
        domain.rho(i,ny-1) = 0 * (Const::eV / (Const::K_b * Const::EV_to_K));
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
                double phi_new = (1/coeff) * ((domain.rho(i,j) - exp(domain.phi(i,j))) * ((domain.L*domain.L)/(domain.LDe*domain.LDe)) 
                + ((domain.phi(i+1,j) + domain.phi(i-1,j))/dx2) + ((domain.phi(i,j+1) + domain.phi(i,j-1))/dy2));
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
                    double R = ((domain.phi(i+1,j) - 2*domain.phi(i,j) + domain.phi(i-1,j))/dx2) + ((domain.phi(i,j+1) - 2*domain.phi(i,j)
                     + domain.phi(i,j-1))/dy2) + (domain.rho(i,j) - exp(domain.phi(i,j))) * ((domain.L*domain.L)/(domain.LDe*domain.LDe));
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
    if(domain.SolverType == "spectral" || domain.bc == "pbc") spectral();
}


///*
/*
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
*/

void FieldSolve::CalculateEfield()
{
    int nx = domain.nx;
    int ny = domain.ny;
    double dx = domain.dx;
    double dy = domain.dy;

    // Electric field along x-direction (E_x = -d(phi)/dx)
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 1; i < nx - 1; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            domain.efx(i,j) = -(domain.phi(i+1,j) - domain.phi(i-1,j)) / (2 * dx);
        }
    }

    if(domain.bc == "open" || domain.bc == "reflective")
    {
        #pragma omp parallel for schedule(static)
        for(int j = 0; j < ny; j++)
        {
            domain.efx(0,j) = -(domain.phi(1,j) - domain.phi(0,j)) / dx;             
            domain.efx(nx-1,j) = -(domain.phi(nx-1,j) - domain.phi(nx-2,j)) / dx; 
        }
    }
    if(domain.bc == "pbc")
    {
        #pragma omp parallel for schedule(static)
        for(int j = 0; j < ny; j++)
        {
            // Since phi(0,j) = phi(nx-1,j), use points just inside the domain for consistency
            domain.efx(0,j) = -(domain.phi(1,j) - domain.phi(nx-2,j)) / (2 * dx);          
            domain.efx(nx-1,j) = domain.efx(0,j); 
        }
    }
    
    // Electric field along y-direction (E_y = -d(phi)/dy)
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 0; i < nx; i++)
    {
        for(int j = 1; j < ny - 1; j++)
        {
            domain.efy(i,j) = -(domain.phi(i,j+1) - domain.phi(i,j-1)) / (2 * dy);
        }
    }
    
    if(domain.bc == "open"  || domain.bc == "reflective")
    {
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < nx; i++)
        {
            domain.efy(i,0) = -(domain.phi(i,1) - domain.phi(i,0)) / dy;                
            domain.efy(i,ny-1) = -(domain.phi(i,ny-1) - domain.phi(i,ny-2)) / dy; 
        }
    }
    if(domain.bc == "pbc")
    {
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < nx; i++)
        {
            // Since phi(i,0) = phi(i,ny-1), use points just inside the domain
            domain.efy(i,0) = -(domain.phi(i,1) - domain.phi(i,ny-2)) / (2 * dy);         
            domain.efy(i,ny-1) = domain.efy(i,0);
        }
    }
}
//*/

/*
void FieldSolve::CalculateEfield()
{
    // Define external electric field components (e.g., uniform field)
    double E_ext_x = 0.0;      // External E-field in x-direction (V/m)
    double E_ext_y = 10;      // External E-field in y-direction (V/m)
  
    for(int i = 1; i < domain.nx - 1; i++)
    {
        for(int j = 0; j < domain.ny; j++)
        {
            domain.efx(i,j) = -(domain.phi(i+1,j) - domain.phi(i-1,j)) / (2 * domain.dx) + E_ext_x;
        }
    }

    if(domain.bc == "open" || domain.bc == "reflective")
    {
        for(int j = 0; j < domain.ny; j++)
        {
            domain.efx(0,j) = -(domain.phi(1,j) - domain.phi(0,j)) / domain.dx + E_ext_x;
            domain.efx(domain.nx-1,j) = -(domain.phi(domain.nx-1,j) - domain.phi(domain.nx-2,j)) / domain.dx + E_ext_x;
        }
    }
    if(domain.bc == "pbc")
    {
        for(int j = 0; j < domain.ny; j++)
        {
            domain.efx(0,j) = -(domain.phi(1,j) - domain.phi(domain.nx-2,j)) / (2 * domain.dx) + E_ext_x;
            domain.efx(domain.nx-1,j) = domain.efx(0,j); // Periodic boundary
        }
    }
    
    for(int i = 0; i < domain.nx; i++)
    {
        for(int j = 1; j < domain.ny - 1; j++)
        {
            domain.efy(i,j) = -(domain.phi(i,j+1) - domain.phi(i,j-1)) / (2 * domain.dy) + E_ext_y;
        }
    }
    
    if(domain.bc == "open" || domain.bc == "reflective")
    {
        for(int i = 0; i < domain.nx; i++)
        {
            domain.efy(i,0) = -(domain.phi(i,1) - domain.phi(i,0)) / domain.dy + E_ext_y;
            domain.efy(i,domain.ny-1) = -(domain.phi(i,domain.ny-1) - domain.phi(i,domain.ny-2)) / domain.dy + E_ext_y;
        }
    }
    if(domain.bc == "pbc")
    {
        for(int i = 0; i < domain.nx; i++)
        {
            domain.efy(i,0) = -(domain.phi(i,1) - domain.phi(i,domain.ny-2)) / (2 * domain.dy) + E_ext_y;
            domain.efy(i,domain.ny-1) = domain.efy(i,0); // Periodic boundary
        }
    }
}
*/
