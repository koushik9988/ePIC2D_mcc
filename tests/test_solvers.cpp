
/*
// tests/test_solvers.cpp
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>
#include <iomanip>
#include "slap.h"

//#ifdef EIGEN_AVAILABLE
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
//#endif

int main()
{

    int nx = 4;
    int ny = 4;
    int n = nx * ny;
    Matrix<double> A(n, n);
    Matrix<double> rhom(nx, ny);

    // Populate matrix A (5-point stencil)
    for (int i = 0; i < n; ++i)
    {
        int row = i / nx;
        int col = i % nx;
        A(i, i) = -4;
        if (col > 0)
            A(i, i - 1) = 1;
        if (col < nx - 1)
            A(i, i + 1) = 1;
        if (row > 0)
            A(i, i - nx) = 1;
        if (row < ny - 1)
            A(i, i + nx) = 1;
    }

    A.display("matrix");

    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            rhom(i, j) = i*j + i+j;
        }
    }

    vec<double> rho(n);
    rho = flatten(rhom);
    vec<double> x0(n); // Initial guess (zero)

    // Solver parameters
    int max_iteration = 100000;
    double tolerance = 1e-5;
    double omega = 1.5; // SOR overrelaxation parameter

    // Store results
    struct Result
    {
        std::string name;
        double time_ms;
        double residual_norm;
    };
    std::vector<Result> results;

    // Helper to compute residual norm
    auto compute_residual = [&](const vec<double>& x) -> double
    {
        vec<double> Ax = A * x;
        vec<double> r(n);
        for (int i = 0; i < n; ++i)
        {
            r(i) = rho(i) - Ax(i);
        }
        return r.norm();
    };

    // Test CG
    {
        auto start = std::chrono::high_resolution_clock::now();
        vec<double> x = cg(A, x0, rho, max_iteration, tolerance);
        auto end = std::chrono::high_resolution_clock::now();
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        results.push_back({"CG", time_ms, compute_residual(x)});
        x.display("cg");
    }

    // Test PCG
    {
        auto start = std::chrono::high_resolution_clock::now();
        vec<double> x = pcg(A, x0, rho, max_iteration, tolerance);
        auto end = std::chrono::high_resolution_clock::now();
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        results.push_back({"PCG", time_ms, compute_residual(x)});
        x.display("pcg");
    }

    // Test Gaussian Elimination
    {
        auto start = std::chrono::high_resolution_clock::now();
        vec<double> x = gausselimination(A, rho);
        auto end = std::chrono::high_resolution_clock::now();
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        results.push_back({"Gaussian Elimination", time_ms, compute_residual(x)});
        x.display("ge");
    }

//#ifdef EIGEN_AVAILABLE
    // Eigen setup
    Eigen::SparseMatrix<double> A_eigen(n, n);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < n; ++i)
    {
        int row = i / nx;
        int col = i % nx;
        triplets.push_back({i, i, -4});
        if (col > 0)
            triplets.push_back({i, i - 1, 1});
        if (col < nx - 1)
            triplets.push_back({i, i + 1, 1});
        if (row > 0)
            triplets.push_back({i, i - nx, 1});
        if (row < ny - 1)
            triplets.push_back({i, i + nx, 1});
    }
    A_eigen.setFromTriplets(triplets.begin(), triplets.end());
    A_eigen.makeCompressed();
    Eigen::VectorXd b_eigen(n);


    std::cout << "\n=== Eigen A_eigen (Full Matrix) ===\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // Access via coeff() returns 0 for non-stored elements
            std::cout << std::setw(5) << std::fixed << std::setprecision(1) << A_eigen.coeff(i, j) << " ";
        }
        std::cout << "\n";
    }

    for (int i = 0; i < n; ++i)
    {
        b_eigen(i) = rho(i); // Use parentheses for consistency
    }
    Eigen::VectorXd x0_eigen = Eigen::VectorXd::Zero(n);

    // Test Eigen CG
    {
        auto start = std::chrono::high_resolution_clock::now();
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
        solver.setMaxIterations(max_iteration);
        solver.setTolerance(tolerance);
        solver.compute(A_eigen);
        Eigen::VectorXd x = solver.solveWithGuess(b_eigen, x0_eigen);
        auto end = std::chrono::high_resolution_clock::now();
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        vec<double> x_slap(n);
        for (int i = 0; i < n; ++i)
        {
            x_slap(i) = x(i);
        }
        results.push_back({"Eigen CG", time_ms, compute_residual(x_slap)});
        x_slap.display("eigen cg");
    }

    // Test Eigen SparseLU
    {
        auto start = std::chrono::high_resolution_clock::now();
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A_eigen);
        Eigen::VectorXd x = solver.solve(b_eigen);
        auto end = std::chrono::high_resolution_clock::now();
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        vec<double> x_slap(n);
        for (int i = 0; i < n; ++i)
        {
            x_slap(i) = x(i);
        }
        results.push_back({"Eigen SparseLU", time_ms, compute_residual(x_slap)});
        x_slap.display("eigen LU");
    }

    // Generate report
    std::cout << "\n=== Solver Performance Report ===\n";
    std::cout << std::setw(20) << "Solver" << std::setw(15) << "Time (ms)" << std::setw(15) << "Residual Norm\n";
    std::cout << std::string(50, '-') << "\n";
    for (const auto& res : results)
    {
        std::cout << std::setw(20) << res.name
                  << std::setw(15) << std::fixed << std::setprecision(2) << res.time_ms
                  << std::setw(15) << std::scientific << std::setprecision(4) << res.residual_norm
                  << "\n";
    }

    return 0;
}
*/

// tests/test_solvers.cpp
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>
#include <iomanip>
#include <random>
#include "slap.h" // Should be "linalg.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

int main()
{
    int nx = 30;
    int ny = 30;
    int n = nx * ny;
    Matrix<double> A(n, n);

    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-10.0, 10.0); // Random values between -10 and 10
    std::uniform_int_distribution<> idx_dis(0, n - 1); // Random indices

    // Populate A with random sparse entries (~5 non-zeros per row)
    const int non_zeros_per_row = 5;
    for (int i = 0; i < n; ++i) {
        std::vector<int> cols;
        A(i, i) = dis(gen); // Ensure diagonal is non-zero
        cols.push_back(i);
        for (int k = 0; k < non_zeros_per_row - 1; ++k) { // Add 4 more non-zeros
            int j = idx_dis(gen);
            if (std::find(cols.begin(), cols.end(), j) == cols.end()) { // Avoid duplicates
                A(i, j) = dis(gen);
                cols.push_back(j);
            }
        }
    }
    //A.display("matrix");

    // Random right-hand side
    vec<double> rho(n);
    for (int i = 0; i < n; ++i) {
        rho(i) = dis(gen);
    }
    vec<double> x0(n); // Zero initial guess

    // Solver parameters
    int max_iteration = 20000;
    double tolerance = 1e-5;

    // Store results
    struct Result {
        std::string name;
        double time_ms;
        double residual_norm;
    };
    std::vector<Result> results;

    // Helper to compute residual norm
    auto compute_residual = [&](const vec<double>& x) -> double {
        vec<double> Ax = A * x;
        vec<double> r(n);
        for (int i = 0; i < n; ++i) {
            r(i) = rho(i) - Ax(i);
        }
        return r.norm();
    };

    // Test CG
    {
        auto start = std::chrono::high_resolution_clock::now();
        vec<double> x = gausselimination(A, rho);
        auto end = std::chrono::high_resolution_clock::now();
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        results.push_back({"CG", time_ms, compute_residual(x)});
        //x.display("cg");
    }

    // Test PCG
    {
        auto start = std::chrono::high_resolution_clock::now();
        vec<double> x = gausselimination(A, rho);
        auto end = std::chrono::high_resolution_clock::now();
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        results.push_back({"PCG", time_ms, compute_residual(x)});
        //x.display("pcg");
    }

    // Test Gaussian Elimination
    {
        auto start = std::chrono::high_resolution_clock::now();
        vec<double> x = gausselimination(A, rho);
        auto end = std::chrono::high_resolution_clock::now();
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        results.push_back({"Gaussian Elimination", time_ms, compute_residual(x)});
        x.display("ge");
    }

    // Eigen setup with identical random matrix
    Eigen::SparseMatrix<double> A_eigen(n, n);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double value = A(i, j);
            if (value != 0) {
                triplets.push_back({i, j, value});
            }
        }
    }
    A_eigen.setFromTriplets(triplets.begin(), triplets.end());
    A_eigen.makeCompressed();


    Eigen::VectorXd b_eigen(n);
    for (int i = 0; i < n; ++i) {
        b_eigen(i) = rho(i);
    }
    Eigen::VectorXd x0_eigen = Eigen::VectorXd::Zero(n);

    auto compute_eigen_residual = [&](const Eigen::VectorXd& x) -> double {
        Eigen::VectorXd Ax = A_eigen * x;
        Eigen::VectorXd r = b_eigen - Ax;
        return r.norm();
    };

    // Test Eigen CG
    {
        auto start = std::chrono::high_resolution_clock::now();
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
        solver.setMaxIterations(max_iteration);
        solver.setTolerance(tolerance);
        solver.compute(A_eigen); // Void return
        if (solver.info() != Eigen::Success) {
            std::cerr << "Eigen CG: Matrix factorization failed\n";
            return 1;
        }
        Eigen::VectorXd x = solver.solveWithGuess(b_eigen, x0_eigen);
        auto end = std::chrono::high_resolution_clock::now();
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        vec<double> x_slap(n);
        for (int i = 0; i < n; ++i) {
            x_slap(i) = x(i);
        }
        double eigen_residual = compute_eigen_residual(x);
        std::cout << "Eigen CG eigen_residual: " << eigen_residual << "\n";
        results.push_back({"Eigen CG", time_ms, compute_residual(x_slap)});
        //x_slap.display("eigen cg");
    }

    // Test Eigen SparseLU
    {
        auto start = std::chrono::high_resolution_clock::now();
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A_eigen); // Void return
        if (solver.info() != Eigen::Success) {
            std::cerr << "Eigen SparseLU: Matrix factorization failed\n";
            return 1;
        }
        Eigen::VectorXd x = solver.solve(b_eigen);
        auto end = std::chrono::high_resolution_clock::now();
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        vec<double> x_slap(n);
        for (int i = 0; i < n; ++i) {
            x_slap(i) = x(i);
        }
        double eigen_residual = compute_eigen_residual(x);
        std::cout << "Eigen SparseLU eigen_residual: " << eigen_residual << "\n";
        results.push_back({"Eigen SparseLU", time_ms, compute_residual(x_slap)});
        x_slap.display("eigen LU");
    }

    std::cout << "\n=== Solver Performance Report ===\n";
    std::cout << std::setw(20) << "Solver" << std::setw(15) << "Time (ms)" << std::setw(15) << "Residual Norm\n";
    std::cout << std::string(50, '-') << "\n";
    for (const auto& res : results) {
        std::cout << std::setw(20) << res.name
                  << std::setw(15) << std::fixed << std::setprecision(2) << res.time_ms
                  << std::setw(15) << std::scientific << std::setprecision(4) << res.residual_norm
                  << "\n";
    }

    return 0;
}