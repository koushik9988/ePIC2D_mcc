#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cstring>
#include <fftw3.h> 
#include "slap.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

void SpectralSolver(Matrix<double> &rho_in, Matrix<double> &phi_out, int nx, int ny, double Lx, double Ly)
{
    int Nh = (ny / 2) + 1;
    double norm = 1.0 / (nx * ny);

    // Allocate real-space and Fourier-space arrays
    double* rho = fftw_alloc_real(nx * ny);
    double* phi = fftw_alloc_real(nx * ny);
    fftw_complex* rho_k = fftw_alloc_complex(nx * Nh);
    fftw_complex* phi_k = fftw_alloc_complex(nx * Nh);

    std::memset(rho, 0, sizeof(double) * nx * ny);
    std::memset(phi, 0, sizeof(double) * nx * ny);
    std::memset(rho_k, 0, sizeof(fftw_complex) * nx * Nh);
    std::memset(phi_k, 0, sizeof(fftw_complex) * nx * Nh);

    // Copy and normalize rho_in to rho array
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            rho[i * ny + j] = -rho_in(i, j);
        }
    }

    // Create FFTW plans
    fftw_plan forward = fftw_plan_dft_r2c_2d(nx, ny, rho, rho_k, FFTW_ESTIMATE);
    fftw_plan backward = fftw_plan_dft_c2r_2d(nx, ny, phi_k, phi, FFTW_ESTIMATE);

    // Forward FFT
    fftw_execute(forward);

    // Solve Poisson equation in Fourier space
    for (int i = 0; i < nx; ++i)
    {
        int kx_index = (i <= nx / 2) ? i : i - nx;
        double kx = 2.0 * M_PI * kx_index / Lx;

        for (int j = 0; j < Nh; ++j)
        {
            double ky = 2.0 * M_PI * j / Ly;
            double k2 = kx * kx + ky * ky;

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

    // Backward FFT
    fftw_execute(backward);

    // Fill phi_out with normalized data
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            phi_out(i, j) = phi[i * ny + j] * norm;
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

int main()
{
    int nx = 1024;
    int ny = 1024;
    double Lx = 1.0;
    double Ly = 1.0;

    Matrix<double> rho(nx, ny);
    Matrix<double> phi(nx, ny);
    Matrix<double> phi_exact(nx, ny);
    Matrix<double> error(nx, ny);

    std::vector<double> xgrid(nx);
    std::vector<double> ygrid(ny);
    for (int i = 0; i < nx; ++i) xgrid[i] = i * Lx / nx;
    for (int j = 0; j < ny; ++j) ygrid[j] = j * Ly / ny;

    // Initialize rho with sin(2pi x) * sin(2pi y)
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            rho(i, j) = std::sin(2.0 * M_PI * xgrid[i] / Lx) * std::sin(2.0 * M_PI * ygrid[j] / Ly);
        }
    }

    // Solve Poisson
    SpectralSolver(rho, phi, nx, ny, Lx, Ly);

    // Compute exact solution
    double coeff = 1.0 / ( (2*M_PI/Lx)*(2*M_PI/Lx) + (2*M_PI/Ly)*(2*M_PI/Ly) );
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            phi_exact(i, j) = coeff * std::sin(2.0 * M_PI * xgrid[i] / Lx) * std::sin(2.0 * M_PI * ygrid[j] / Ly);
            error(i, j) = phi(i, j) - phi_exact(i, j);
        }
    }

    vec<double> error_flatten = flatten(error);
    std::cout<<error_flatten.norm()<<std::endl;

    // Prepare data for plotting
    std::vector<std::vector<double>> phi_data(ny, std::vector<double>(nx));
    std::vector<std::vector<double>> phi_exact_data(ny, std::vector<double>(nx));
    std::vector<std::vector<double>> error_data(ny, std::vector<double>(nx));

    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            phi_data[j][i] = phi(i,j);
            phi_exact_data[j][i] = phi_exact(i,j);
            error_data[j][i] = error(i,j);
        }
    }

    // Plot numerical solution
    plt::figure(1);
    plt::imshow(phi_data, "coolwarm", "lower", "bilinear");
    //plt::colorbar();
    plt::title("Numerical phi");
    plt::xlabel("x");
    plt::ylabel("y");

    // Plot exact solution
    plt::figure(2);
    plt::imshow(phi_exact_data, "coolwarm", "lower", "bilinear");
    //plt::colorbar();
    plt::title("Exact phi");
    plt::xlabel("x");
    plt::ylabel("y");

    // Plot error
    plt::figure(3);
    plt::imshow(error_data, "coolwarm", "lower", "bilinear");
    //plt::colorbar();
    plt::title("Error (Numerical - Exact)");
    plt::xlabel("x");
    plt::ylabel("y");

    plt::show();
}
