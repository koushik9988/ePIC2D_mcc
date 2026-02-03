# Electrostatic 2D-3V Particle-in-Cell (PIC) Code

This repository contains an electrostatic 2D-3V Particle-in-Cell (PIC) code. The code is in development stage and README incomplete.




## Requirements
- Python3 : Required for data processing, and data visualization. 
- python3-dev : Provides Python development headers needed for matplotlibcpp.
- GNU C++ compiler / clang
- [CMake](https://cmake.org/)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [matplotlibcpp](https://github.com/lava/matplotlib-cpp)
- [Git](https://git-scm.com/)
- Matplotlib
- NumPy
- Scipy


### Installation
1. Clone the repository:
    ```bash
    git clone https://github.com/koushik9988/ePIC2D_mcc.git
    ```

2. Navigate to the directory:
    ```bash
    cd ePIC_mcc
    ```

3. Build the code using cmake:
    ```bash
    mkdir build && cd build
    ```
    ```bash
    cmake ..
    ```
    ```bash
    cmake --build .
    ```

### Running the Code
1. Configure the simulation parameters in the `input.ini` file.
2. Run the code:
The executble will be located in the build directory after building with cmake.
    ```bash
    ./ePIC++ ../inputfiles/input.ini
    ```

# Explanation of `input.ini` File Parameters

The `input.ini` file contains parameters for configuring the simulation. Each section corresponds to a different aspect of the simulation setup.

## `[file]`

| Parameter | Description |
|----------|-------------|
| `output` | Path to directory where simulation output data in hdf5 format is saved |

## `[time]`
| Parameter | Description |
|----------|-------------|
| `NUM_TS` | Total number of time steps to run |
| `DT_coeff` | Time step as a fraction of (1/frequency): `dt = DT_coeff / ω_pe` or `dt = DT_coeff / ω_pi` |

## `[diagnostics]`
| Parameter | Description |
|----------|-------------|
| `write_interval` | Interval (in steps) to write density/field data |
| `write_interval_phase` | Interval to write phase-space data |
| `write_diagnostics` | Interval to output diagnostics in screen (energies, phase plot etc.) |
| `write_flag` | What data to write in disk: `0 = none`, `1 = all`, `2 = fields only`, `3 = phase-space only` |
| `save_fig` | Save plots as images (`1 = yes`, `0 = no`) (deprecated) |
| `sub_cycle_interval` | Frequency of ion sub-cycling  |
| `precision` | flaoting point precision level |
| `diagtype` | Type of diagnostics :`off(print just time step)`, `basic (print just time and max_phi)` ,`full (print evryting with live plot)`|

## `[PlotFlags]`

| Flag | Description |
|------|-------------|
| `phase_space` | Phase-space plot: `1 = x-vx`, `2 = y-vy` |
| `config_space` | Configuration plot: `1 = x-y`, `2 = vx-vy` ( `vx-vy` isnot implemented) |
| `electric_field` | Electric field plot: `1 = Ex`, `2 = Ey` |
| `potential_field` | Enable potential field plotting |
| `density_contour` | Plot species density contour |
| `ke_components` | Plot kinetic energy components (x, y, z) |
| `total_energy` | Plot total energy evolution |
| `species_index` | Index of species to use for phase-space and density plots starting from index 0 as in species section of the input file (e.g 0 = electrons, 1 = ion etc) |
 
## `[domain]`

| Parameter | Description |
|----------|-------------|
| `nx` | Number of grid points along x axis |
| `nx` | Number of grid points along y axis |
| `x0` | Origin coordinate of the domain |
| `y0` | Origin coordinate of the domain |

## `[Grids]`
| Property         | Description |
|------------------|-------------|
|`grid_count`      |number of grids (`first grid paramter start with grid_0 , second with grid_1 and so on`)|
| `type`           | Grid type. Supported: `reactconduct` (biased), `circular` (biased), (will implement more later)|
| `min_x`, `max_x` | X-range (cell indices) of the reactconduct grid |
| `min_y`, `max_y` | Y-range (cell indices) of reactconduct the grid |
| `grid_voltage`   | Voltage applied to the grid |
| `x_center` , `y_center`  | center cocordinate of circular grid|
| `electrode_number`  | number of electrode (circular shape) the circular grid consist of|
| `grid_radius`  | radious of circular grid|
| `electrode_radius`  | radious of electrode|

```
grid_count = 3
grid_0.type = reactconduct
grid_0.min_x = 0
grid_0.min_y = 30
grid_0.max_x = 22
grid_0.max_y = 31
grid_0.grid_voltage = 100 # in volt

grid_1.type = reactconduct
grid_1.min_x = 0
grid_1.min_y = 30
grid_1.max_x = 22
grid_1.max_y = 31
grid_1.grid_voltage = 100 # in volt

grid_2.type = circular
grid_2.x_center = 40   
grid_2.y_center = 25
grid_2.electrode_number = 10 
grid_2.grid_radius = 8        
grid_2.electrod_radius = 1   
grid_2.electrode_voltage = -100  
```
## `[Emitters]`
| Property            | Description |
|---------------------|-------------|
| `x0`, `y0`          | Starting coordinate of emitter line (in grid units) |
| `x1`, `y1`          | Ending coordinate of emitter line (in grid units) |
| `temp`              | Emission temperature (not implemented if `0`) |
| `numparticle`       | Number of particles injected per time step |
| `vdx`, `vdy`        | Initial velocity components for injected particles |
| `species_idx1`      | Index of the first species to emit |
| `species_idx2`      | Index of the 2nd species to emit|
|if both index `species_idx1` and `species_idx2` are same then emitter will emit only emit that species (no double emission)|
```
count = 0  

emitter_0.x0 = 1
emitter_0.y0 = 20
emitter_0.x1 = 1
emitter_0.y1 = 30
emitter_0.temp = 0
emitter_0.numparticle = 1
emitter_0.vdx = 2.0
emitter_0.vdy = 0.0
emitter_0.species_idx1 = 0
emitter_0.species_idx2 = 1
```

## `[normalization]`

| Parameter | Description |
|-----------|-------------|
| `norm_scheme` | Type of normalization (`1 = electron scale`, `2 = ion scale`, `3 = subcycling`, `4 = mixed scale`, `5 = mixed scale`) |
| `vel_norm_scheme` | Velocity normalization type (`1 = vthe`, `2 = vthi`, `3 = vcs(ion-acoustic speed)`) |
| `lenght_scale` | User defined Characteristic length scale |
| `time_scale` | User Defined Characteristic time scale (`1/ω_pe`) |
| `energy_scale` | User defined Characteristic energy scale |

## `[simulation]`

| Parameter | Description |
|-----------|-------------|
| `shapefunction` | Interpolation scheme: `NGP` and  `CIC` |
| `push_parallal` | Use parallel particle pushing (`1 = enabled`) |
| `deposit_parallal` | Use parallel charge deposition (`1 = enabled`) |
| `density` | Plasma density |
| `bc` | Boundary condition: `pbc` (periodic) or `open` |
| `see_rate` | Secondary electron emission rate *(not implemented)* |
| `tempwall` | Wall temperature *(not implemented)* |
| `ionfixed` | Fixed background ions (`1 = yes`, `0 = no`) |

## `[solver]`

| Parameter | Description |
|----------|-------------|
| `solvertype` | Solver method: `gs`, `pcg`,`spectral` `(use spectral for periodic boundary)`|
| `tolerance` | Convergence tolerance (for iterative solvers) |
| `max_iteration` | Maximum iterations (for iterative solvers) |

## `[collision]`
| Parameter | Description |
|----------|-------------|
| `elastic` | Enable elastic collisions (`true` and `false`) |
| `excitation` | Enable excitation collisions |
| `ionization` | Enable ionization |
| `GAS_DENSITY` | Neutral gas density (e.g., `1e20`) |
| `collgroup ` | particle collision group in a pair of two |
```
collgroup = 12,21 means that first species will collide with a gas which is neutral form of 2nd species (whichis generallu an ion)
and third species will collide with same gas which is neutral form of 2nd species 
```

## `[species]`

Each line represents a species and its properties in the following format:

  ```
  name, mass, number_of_particles, temperature, charge_sign, normalized density (w.r.t electron density), streaming_velocity along x direction, streaming_velocity along y direction, load_type
  ```
  
Example species configuration:
  
  ```
  electron, 9.10938215E-31, 50000, 1, -1, 1, -10,0, uniform
  ion, 6.63352090e-26, 50000, 0, 1, 0, 0, 0,uniform
  beam, 9.10938215E-31, 50000, 1, -1, 1, 10,0, uniform
  ```
(Note : Electron should be in the first line and Ion should be in the 2nd line and all other species will go after that.)


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
  



