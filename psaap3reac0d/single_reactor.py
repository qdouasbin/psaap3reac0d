import os
import pprint

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import toml


def single_reactor_simulation(params):
    """
    Single reactor network for a rocket combustion chamber
    """

    print_params(params)

    # Get the simulation parameters
    _mech = params['Case']['chemical_mechanism']
    _fuel = params['FluidConditions']['Fuel']
    _oxidizer = params['FluidConditions']['Oxidizer']
    _ignition = params['Ignition']
    _exhaust = params['FluidConditions']['Exhaust']
    _geom_ch = params["Geometry"]['CombustionChamber']
    _geom_nozzle = params['Geometry']['ExhaustNozzle']
    _solver = params['Solver']
    _out_dir = params['IO']['output_dir']

    if params['InitialSolution']['type'] == 'Exhaust':
        _initial_conditions = params['FluidConditions']['Exhaust']
    else:
        raise NotImplementedError(
            "The only initial solution implemented is for a combustor "
            "initially filled with the temperature, pressure and composition of the surroundings.")

    # Additional variables derived from input data
    _geom_ch['area'] = np.pi * _geom_ch['radius'] ** 2.
    _geom_nozzle['area'] = np.pi * _geom_nozzle['radius'] ** 2.

    # Fuel
    fuel_gas = ct.Solution(_mech)
    fuel_gas.TPX = _fuel['temperature'], _fuel['pressure'], _fuel['molar_frac']
    reservoir_fuel = ct.Reservoir(fuel_gas)

    # Oxidizer
    ox_gas = ct.Solution(_mech)
    ox_gas.TPX = _oxidizer['temperature'], _oxidizer['pressure'], _oxidizer['molar_frac']
    reservoir_oxidizer = ct.Reservoir(ox_gas)

    # Igniter --> H+ ions
    igniter_gas = ct.Solution(_mech)
    igniter_gas.TPX = _ignition['temperature'], _ignition['pressure'], _ignition['molar_frac']
    igniter = ct.Reservoir(igniter_gas)

    # Create the combustor, and fill it initial conditions
    combustor_gas = ct.Solution(_mech)
    combustor_gas.TPX = _initial_conditions['temperature'], \
                        _initial_conditions['pressure'], \
                        _initial_conditions['molar_frac']
    combustor = ct.IdealGasReactor(combustor_gas)
    combustor.volume = _geom_ch['length'] * _geom_ch['area']

    # Create a reservoir for the exhaust
    exhaust_gas = ct.Solution('gri30.xml')
    exhaust_gas.TPX = _exhaust['temperature'], \
                      _exhaust['pressure'], \
                      _exhaust['molar_frac']
    exhaust = ct.Reservoir(exhaust_gas)

    # Create mass flow controller (Fuel reservoir --> combustion chamber)
    mfc_fuel = ct.MassFlowController(reservoir_fuel, combustor, mdot=_fuel['mdot'])

    # Create mass flow controller (Oxidizer reservoir --> combustion chamber)
    mfc_oxidizer = ct.MassFlowController(reservoir_oxidizer, combustor, mdot=_oxidizer['mdot'])

    # Create mass flow controller (igniter reservoir --> combustion chamber)
    # mdot = mdot(t)
    def gaussian_pulse_H_igniter_mdot(t):
        """
        Gaussian pulse of H+ ions centered around t0
        :param t: time of the simulation
        :return: mdot of H+ [kg/s]
        """
        amplitude = _ignition['maximum_mass_flux']
        t0 = _ignition['time']
        fwhm = _ignition['fwhm']
        return amplitude * np.exp(-(t - t0) ** 2 * 4 * np.log(2) / fwhm ** 2)

    mfc_igniter = ct.MassFlowController(igniter, combustor, mdot=gaussian_pulse_H_igniter_mdot)

    def nozzle_nasa_mdot(t):
        """
        Compute the efflux of the combustor via simple isentropic relations.
        The efficiency of the nozzle is modeled using a discharge coefficient
        :param t: time
        :return: mass flux at the outlet (efflux) of the combustion chamber.
        """
        # Avoid a unity pressure ratio at startup
        epsilon = 1.0
        _rho = combustor.thermo.density_mass
        _press = combustor.thermo.P
        _area = _geom_nozzle['area']
        _gamma_s = combustor.thermo.cp_mass / combustor.thermo.cv_mass
        _temp = combustor.T
        _r_gas_specific = combustor.thermo.cp_mass - combustor.thermo.cv_mass
        _p_ratio = _exhaust['pressure'] / (combustor.kinetics.P + epsilon)

        power_1 = 2.0 / _gamma_s
        power_2 = (_gamma_s + 1.0) / (_gamma_s)
        pressure_term = _p_ratio ** power_1 - _p_ratio ** power_2
        _sqrt_term = 2. * _gamma_s * _r_gas_specific * _temp / (_gamma_s - 1.0)
        _sqrt_term *= pressure_term
        assert (_sqrt_term >= 0.0)
        _sqrt_term = np.sqrt(_sqrt_term)

        _mdot_unchoked = _rho * _area * _sqrt_term

        power = (_gamma_s + 1.0) / (_gamma_s - 1.0)
        _gamma_term = (2. / (_gamma_s + 1.0)) ** power
        _sqrt = np.sqrt(_gamma_s * _r_gas_specific * _temp * _gamma_term)
        _mdot_choked = _rho * _area * _sqrt

        _mdot = 0.0
        _p_crit_downstream = _press * (2. / (_gamma_s + 1.0)) ** (_gamma_s / (_gamma_s - 1.0))

        if _p_crit_downstream < _exhaust['pressure']:
            # Unchoked conditions
            _mdot = _mdot_unchoked
        else:
            # Choked conditions
            _mdot = _mdot_choked

        return _mdot * _geom_nozzle['C_d']

    # Create mass flow controller for the combustor efflux (combustion chamber --> exhaust)
    # mdot = mdot(p_chamber, rho_chamber, p_exhaust, area_nozzle_throat, C_d)
    mfc_nozzle = ct.MassFlowController(combustor, exhaust, mdot=nozzle_nasa_mdot)

    # Create a single reactor network
    sim = ct.ReactorNet([combustor])

    # Solve and store solution
    t_now = _solver['t_init']
    t_final = _solver['t_end']
    temp_prev = combustor.T
    t_prev = t_now

    # Store results in a dictionary
    res = {}
    res['time'] = []
    res['T'] = []
    res['P'] = []
    res['rho'] = []
    res['residence_time'] = []
    res['mdot_out'] = []
    res['mdot_fuel'] = []
    res['mdot_ox'] = []
    res['mdot_igniter'] = []
    res['mass'] = []

    while t_now < t_final:
        t_now = sim.step()
        temp_now = combustor.T

        # Store data point at least every 5 Kelvin or every millisecond
        if abs(temp_now - temp_prev) > 5 or t_now - t_prev > 1e-3:
            t_prev = t_now
            temp_prev = temp_now
            t_res = combustor.mass / mfc_nozzle.mdot(t_now)
            cur_state = ct.Solution('gri30.xml')
            _cur_state_rho = combustor.mass / combustor.volume
            _cur_state_T = combustor.T
            _cur_state_P = combustor.kinetics.P
            print(' > t = %e [s]    T = %e [K]    p = %e [Pa]' % (t_now, temp_now, _cur_state_P))

            res['time'].append(t_now)
            res['T'].append(temp_now)
            res['P'].append(_cur_state_P)
            res['rho'].append(_cur_state_rho)
            res['residence_time'].append(t_res)
            res['mdot_out'].append(nozzle_nasa_mdot(t_now))
            res['mdot_fuel'].append(_fuel['mdot'])
            res['mdot_ox'].append(_oxidizer['mdot'])
            res['mdot_igniter'].append(gaussian_pulse_H_igniter_mdot(t_now))
            res['mass'].append(combustor.mass)

    df = pd.DataFrame.from_dict(res)
    df['mdot_in'] = df['mdot_fuel'] + df['mdot_ox'] + df['mdot_igniter']

    # Create output folder if it does not exists
    if not os.path.exists(_out_dir):
        os.makedirs(_out_dir)

    # Output data as a csv file
    df.to_csv(os.path.join(_out_dir, '%s_%s.csv' % (params['Case']['name'], _ignition['type'])), index=False)

    return df


def print_params(params):
    """Print parameters"""
    print("\n > -------- Parameters -------- ")
    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(params)
    print(" > ---------------------------- \n")


if __name__ == "__main__":
    # load TOML input file as a dictionary
    input_file = '../tests/BKD_LP4.toml'
    params = toml.load(input_file)

    print(params)

    result = single_reactor_simulation(params)
