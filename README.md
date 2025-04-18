# Electrostatic 1D Particle-in-Cell (PIC) Code

This repository contains an electrostatic 2D-3V Particle-in-Cell (PIC) code. The code is in development stage and README incomplete.




## Requirements
- Python3
- GNU C++ compiler(g++)
- GNU make
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- git
- Matplotlib
- NumPy
- Scipy


### Installation
1. Clone the repository:
    ```bash
    git clone https://github.com/koushik9988/particle-in-cell.git
    ```

2. Navigate to the directory:
    ```bash
    cd your_repository
    ```

3. Build the code using make:
    ```bash
    make clean
    ```
    ```bash
    make all
    ```

### Running the Code
1. Configure the simulation parameters in the `input.ini` file.
2. Run the code:
    ```bash
    g++ ./pic
    ```

# Explanation of input.ini file Parameters

The input file `input.ini` contains parameters for configuring the simulation. Each section corresponds to a different aspect of the simulation setup.

## `[file]`

- **output**: Specifies the directory where the output data will be stored.

## `[time]`

- **NUM_TS**: Total number of time steps for the simulation.

## `[diagnostics]`

- **write_interval**: Interval for writing density and field data in result.txt file.
- **write_interval_phase**: Interval for writing phase-space data file.
- **write_diagnostics**: Interval for writing diagnostic outputs.
- **DT_coeff**: Time step.
  
 $$dt = DT_{\text{coeff}} \frac{1}{\omega_{\text{pe}}}$$

 Normalized time step

 $$dt = DT_{\text{coeff}} \frac{1}{\omega_{\text{pe}}}{\omega_{pe}}$$

  

- **write_flag**: Flag for controlling data writing: 
  - 1: Write both phase and field data.
  - 2: Write only field data.
  - 3: Write only phase data.
  - 0: Write no data.

## `[domain]`

- **NC**: Number of cells in the domain.
- **x0**: Initial position of the domain.

## `[population]`

- **nParticlesE**: Number of electrons loaded into the domain.
- **nParticlesI**: Number of ions loaded into the domain.
- **nParticlesN**: Number of negative ions loaded into the domain.
- **nParticlesB**: Number of background particles.
- **tempE**: Temperature of electrons.
- **tempI**: Temperature of ions.
- **tempN**: Temperature of negative ions.
- **tempB**: Temperature of background particles.
- **massE**: Mass of electrons.
- **massI**: Mass of ions.
- **massN**: Mass of negative ions.
- **massB**: Mass of background particles.

## `[simulation]`
- **number_of_species**: Number of species.
- **v_i**: Ion streaming velocity.
- **v_e**: Electron streaming velocity.
- **v_b**: Beam particle streaming velocity.
- **v_n**: Negative ion streaming velocity.
- **density**: Plasma density.
- **alpha**: Fraction of negative ion to background positive ion
- **beta**: fraction of negative beam to background negative ion 
- **bc**: Boundary condition.
   - 1: pbc for periodic boundary.
   - 2: open for open boundary condition.

 # Data processing and visualization
 1. Plot kinetic enegy ,potential enegy and total enegy
     ```bash
    python3 ke_plot.py ../name_of_outputfolder
    ```
 2. Plot dispersion
     ```bash
    python3 dispersion.py ../name_of_outputfolder
    ```
 3. Plot/Animate phase-space and potential data
     ```bash
    python3 phase_pot_plot.py ../name_of_outputfolder
    ```

## Contributors
- Rakesh Moulick
- Kaushik Kalita
  



