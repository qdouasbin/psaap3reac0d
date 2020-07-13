# Reduced Order Modeling for the transient ignition of combustion chambers

We want to monitor the pressure and temperature traces of a rocket combustor during (laser) ignition.
To do so, we use a Cantera network built from simple 0D components (reservoirs, mass flow controllers, 0D reactors, ...).

## Theory

A theory guide is available ![here](https://cantera.org/science/reactors.html).

## Installation

 1. Cantera installation using Conda: ![https://cantera.org/install/conda-install.html](https://cantera.org/install/conda-install.html)
 2. git clone of the repository

Todo: Cantera install + git clone + python setup.py install/devel

## How to use it

The input files are in the TOML format. 
This allows to handle the parameters of the solver as a nested dictionary, which provides an easy API 
to run ensemble runs or parameteric studies.


Here is an example:

```TOML
# Input file for the BKD model rocket (DLR).
# Here, we study the operating condition called LP4.
# See this paper for more details: https://www.sciencedirect.com/science/article/pii/S1540748916301006 

[Case]
name = 'BKD_LP4'
# chemical mechanism used for all components
chemical_mechanism = 'gri30.xml'

[Geometry]
    [Geometry.CombustionChamber]
    # [m]
    length = 20e-2
    # [m]
    radius = 4e-2

    [Geometry.ExhaustNozzle]
    # radius at throat [m]
    radius = 0.025
    # Discharge coefficient (e.g. the efficiency of the nozzle) [-]. Typically, C_d belongs to [0.95, 1]
    # Note, Discharge coefficient taken from this report: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19750006921.pdf
    C_d = 0.99

[FluidConditions]
    [FluidConditions.Fuel]
    # Temperature [K]
    temperature = 96.0
    # molar mass fraction string (Cantera format, e.g. 'CH4:0.05,H2:0.001'),
    molar_frac = 'H2:1.0'
    # pressure in fuel dome [Pa]
    pressure = 103e5
    # Mass flux (influx) [kg/s]
    mdot = 0.96

    [FluidConditions.Oxidizer]
    # Temperature [K]
    temperature = 111.0
    # molar mass fraction string (Cantera format, e.g. 'CH4:0.05,H2:0.001'),
    molar_frac = 'O2:1.0'
    # pressure in oxidizer dome [Pa]
    pressure = 94e5
    # Mass flux (influx) [kg/s]
    mdot = 5.75

[FluidConditions.Exhaust]
# Temperature [K]
temperature = 300.0
# molar mass fraction string (Cantera format, e.g. 'CH4:0.05,H2:0.001'),
molar_frac = 'N2:1.0'
# pressure in oxidizer dome [Pa]
pressure = 94e5

[Ignition]
type = 'Gaussian_pulse_H+'
temperature = 100.0
pressure = 100.0e5
molar_frac = 'H:1.0'
# Full Width at half maximum
fwhm = 1e-3
# Time at the center of the Gaussian pulse [s]
time = 2e-2
maximum_mass_flux = 0.5

[InitialSolution]
type = 'Exhaust'

[Solver]
# Time at the begining of the simulation [s].
t_init = 0.0
# End time of the simulation, [s]
t_end = 3e-2

[IO]
# Output directory (right results as csv file)
output_dir = 'results'
```

