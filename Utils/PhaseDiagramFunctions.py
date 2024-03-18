"""
Title: Spacetime magnetic hopfions: from internal excitations and braiding of skyrmions

Authors and Affiliations: Ross Knapman[1,2,3], Timon Tausendpfund[1], Sebastián A. Díaz[2,4], Karin Everschor-Sitte[2,3]
    1 Institute of Physics, Johannes Gutenberg University Mainz, 55128 Mainz, Germany
    2 Faculty of Physics, University of Duisburg-Essen, 47057 Duisburg, Germany
    3 Center for Nanointegration Duisburg-Essen (CENIDE), University of Duisburg-Essen, 47057 Duisburg, Germany
    4 Department of Physics, University of Konstanz, 78457 Konstanz, Germany

Purpose: Collection of functions used to create the phase diagrams shown in the manuscript.

Written 2023 by Ross Knapman <ross.knapman@uni-due.de>
"""

import re
import os
import numpy as np
import matplotlib
from Utils.Fits import *


def get_E0_and_omega_arrays(R_directory):

    """
    For a given directory containing the directory of evolution arrays of R
    in the phase diagram, obtain the amplitude and frequency arrays.

    Args:
        R_directory (str): The directory in which the integrated collective coordinate integration of
            skyrmion radius R are stored
    
    Returns:
        The values of E0 and omega for which the integrations were carried out by having read the filenames
        in this directory.

    """

    # Obtain list of [[E0_1, omega_1], [E0_2, omega_2], ...] for all files in directory
    extracted_nos = list(map(lambda x: re.findall(r'[0-9]+\.[0-9]+', x), os.listdir(R_directory)))

    # Extract only the E0 values from the above list
    E0_list = [numbers[0] for numbers in extracted_nos]
    
    # Discard duplicate values from above E0 list, sort it, convert to NumPy array
    E0_values = np.array(sorted(list(set(list(E0_list)))), dtype=float)
    
    # Make sure the data is equally spaced as it should be
    assert np.all(np.isclose(np.diff(E0_values), np.diff(E0_values)[0]))
    
    # Same as with E0 above but with omega
    omega_list = [numbers[1] for numbers in extracted_nos]
    omega_values = np.array(sorted(list(set(list(omega_list)))), dtype=float)
    assert np.all(np.isclose(np.diff(omega_values), np.diff(omega_values)[0]))

    return E0_values, omega_values


def get_coordinates(E0, omega, cc_directory):

    """
    For a given amplitude and omega value in the phase diagram, for a given calculation
    (specified by cc_directory), obtain the collective coordinate arrays.

    Args:
        E0 (float): The electric field amplitude
        omega (float): The electric field oscillation angular frequency
        cc_directory (str): The directory in which the collective coordinate data is contained

    Returns:
        The NumPy arrays of skyrmion radius and helicity loaded from the files.

    """

    filename = "E0{:.2f}".format(E0) + "omega{:.2f}".format(omega) + ".npy"
    R = np.load(cc_directory + '/R/' + filename)
    eta = np.load(cc_directory + '/eta/' + filename)
    return R, eta


def get_winding_numbers(R, eta, times, omega):
    
    """
    Bin coordinate arrays into single periods of Ez and calculate array of winding numbers.

    Args:
        R (ndarray): Array of skyrmion radii from collective coordinate integration
        eta (ndarray): Array of skyrmion helicities from collective coordinate integration
        times (ndarray): Array of times of the collective coordinate integration
        omega (float): Electric field oscillation angular frequency

    Returns:
        Array of winding numbers for each electric field oscillation cycle.
    
    """

    bin_width = 2*np.pi / omega
    bin_width_indices = int(np.round(bin_width / (times[1] - times[0])))

    # Number of cycles being considered (from starting point to the end of the simulation, ignore the fraction of a period at the end)
    no_cycles = int(np.floor(len(R) / bin_width_indices))
    winding_nos = np.zeros(no_cycles, dtype=float)

    for i in range(no_cycles):
        start_idx = i*bin_width_indices
        end_idx = (i+1)*bin_width_indices
        winding_nos[i] = np.sum(np.diff(eta[start_idx:end_idx]) / (2*np.pi))

    return winding_nos


def get_skyrmion_energy(cc_directory):

    """
    Get the minimal energy of the skyrmion in the absence of an electric field.

    Args:
        cc_directory (str): The directory in which the collective coordinate data is contained

    Returns:
        The energy for a skyrmion with a radius that minimises its energy.

    """

    Bz, dw_width, E_exchangeFitParams, E_magnetic_integralFitParams, E_electric_integralFitParams, \
            alpha, F_RexFitParams, Gamma11FitParams, Gamma22FitParams, G12FitParams, Xi1FitParams, \
            Xi2FitParams, Xi3FitParams, Xi4FitParams, optimalSkyrmionRadius = read_parameters(cc_directory)

    R_values = np.linspace(0.5, 2.0, 1000)
    E_values = np.zeros_like(R_values, dtype=float)

    for i in range(len(R_values)):
        E_values[i] = total_energy(R_values[i], np.pi, 1., 0., E_exchangeFitParams, E_magnetic_integralFitParams, E_electric_integralFitParams)

    return np.min(E_values)


def total_energy(R, eta, Bz, Ez, E_exchangeFitParams, E_magnetic_integralFitParams, E_electric_integralFitParams):

    """
    Get array of total energies (i.e. exchange + Zeeman + electric contributions) over the course of the collective coordinate integration.

    Args:
        R (ndarray): Array of radii obtained by the collective coordinate integration
        eta (ndarray): Array of helicities obtained by the collective coordinate integration
        Bz (ndarray): Array of magnetic field values during the collective coordinate integration
        Ez (ndarray): Array of electric field values during the collective coordinate integration
        E_exchangeFitParams: Array of fitting parameters for the fit of exchange energy over radius (calculated in CCIntegration.sage and stored in CCParams.txt)
        E_magnetic_integralFitParams: Array of fitting parameters for the fit of Zeeman energy integral over radius (calculated in CCIntegration.sage and stored in CCParams.txt)
        E_electric_integralFitParams: Array of fitting parameters for the fit of electric field energy integral over radius (calculated in CCIntegration.sage and stored in CCParams.txt)

    Returns:
        Array of total energy during the collective coordinate integration

    """

    exchange_energy = inverse_quadratic_linear_fit(R, E_exchangeFitParams[0], E_exchangeFitParams[1], E_exchangeFitParams[2], E_exchangeFitParams[3])
    magnetic_energy = -Bz * quadratic_fit(R, E_magnetic_integralFitParams[0], E_magnetic_integralFitParams[1], E_magnetic_integralFitParams[2])
    electric_energy = -Ez * np.cos(eta) * linear_fit(R, E_electric_integralFitParams[0], E_electric_integralFitParams[1])

    return exchange_energy + magnetic_energy + electric_energy


def get_mean_energy(cc_directory, E0_values, omega_values, Eidx, omegaidx, energy_function, cutoff=0.8):

    """
    Calculate the mean energy for several periods of the electric field oscillation for the last 20%
    of the simulation (this percentage is determined by the cutoff argument).

    Args:
        cc_directory (str): The directory in which the collective coordinate data is contained
        E0_values (ndarray): Array of electric field amplitudes for which the collective coordinate were time-integrated
        omega_values (ndarray): Array of electric field angular frequencies for which the collective coordinate were time-integrated
        Eidx (ndarray): The index within the array E0_values for which the mean energy is to be calculated
        omegaidx (ndarray): The index within the array omega_values for which the mean energy is to be calculated
        energy_function (function): The function outputting the energies (e.g. total_energy)
        cutoff (float): The fraction of the simulation to be excluded from the calculation of the mean energy (e.g. to exclude initial transient)

    Returns:
        Array of mean energies, with a length that depends on the size of the period for that specific omega

    """

    if omegaidx != 0:

        # Electric field amplitude and frequency
        E0 = E0_values[Eidx]
        omega = omega_values[omegaidx]

        times = np.load(cc_directory + '/General/times.npy')
        R_array, eta_array = get_coordinates(E0, omega, cc_directory)

        Bz, dw_width, E_exchangeFitParams, E_magnetic_integralFitParams, E_electric_integralFitParams, \
            alpha, F_RexFitParams, Gamma11FitParams, Gamma22FitParams, G12FitParams, Xi1FitParams, \
            Xi2FitParams, Xi3FitParams, Xi4FitParams, optimalSkyrmionRadius = read_parameters(cc_directory)

        # Array of magnetic and electric fields
        Bz_array = np.full_like(times, Bz)
        Ez_array = E0 * np.cos(omega*times)

        Rcut = R_array[int(cutoff*len(R_array)):]
        etacut = eta_array[int(cutoff*len(eta_array)):]
        Ez_arraycut = Ez_array[int(cutoff*len(times)):]
        Bz_arraycut = Bz_array[int(cutoff*len(times)):]

        bin_width = 2*np.pi / omega
        bin_width_indices = int(np.round(bin_width / (times[1] - times[0])))
        no_cycles = int(np.floor(len(Rcut) / bin_width_indices))

        mean_energies = np.zeros(no_cycles, dtype=float)

        for cycle_idx in range(no_cycles):
            cycle_start_idx = cycle_idx*bin_width_indices
            cycle_end_idx = (cycle_idx+1)*bin_width_indices
            
            to_average = energy_function(Rcut[cycle_start_idx:cycle_end_idx], etacut[cycle_start_idx:cycle_end_idx],
            Bz_arraycut[cycle_start_idx:cycle_end_idx], Ez_arraycut[cycle_start_idx:cycle_end_idx],
            E_exchangeFitParams, E_magnetic_integralFitParams, E_electric_integralFitParams)

            mean_energies[cycle_idx] = np.mean(to_average)

        return np.mean(mean_energies)

    else:
        return 0


def max_radius(cc_directory, E0_values, omega_values, Eidx, omegaidx):

    """
    Obtain the maximum skyrmion radius attained during the collective coordinate evolution.

    Args:
        cc_directory (str): The directory in which the collective coordinate data is contained
        E0_values (ndarray): Array of electric field amplitudes for which the collective coordinate were time-integrated
        omega_values (ndarray): Array of electric field angular frequencies for which the collective coordinate were time-integrated
        Eidx (ndarray): The index within the array E0_values for which the mean energy is to be calculated
        omegaidx (ndarray): The index within the array omega_values for which the mean energy is to be calculated

    Returns:
        The maximum radius attaned during the second half of the collective coordinate time integration (only second half to cut out initial transient)

    """

    R, eta = get_coordinates(E0_values[Eidx], omega_values[omegaidx], cc_directory)

    return np.max(R[len(R)//2:])


def mean_winding_number(cc_directory, E0_values, omega_values, Eidx, omegaidx):

    """
    Obtain the mean number of times the skyrmion's helicity rotates per oscillation period of the electric field.

    Args:
        cc_directory (str): The directory in which the collective coordinate data is contained
        E0_values (ndarray): Array of electric field amplitudes for which the collective coordinate were time-integrated
        omega_values (ndarray): Array of electric field angular frequencies for which the collective coordinate were time-integrated
        Eidx (ndarray): The index within the array E0_values for which the mean energy is to be calculated
        omegaidx (ndarray): The index within the array omega_values for which the mean energy is to be calculated

    Returns:
        The mean number of times the helicity rotates during an electric field oscillation, during the second half of
        the collective coordinate time integration (positive for anticlockwise rotation, negative for clockwise)

    """

    if omegaidx != 0:  # Otherwise get zero division error when calculation time period

        R, eta = get_coordinates(E0_values[Eidx], omega_values[omegaidx], cc_directory)
        times = np.load(cc_directory + '/General/times.npy')

        Rcut = R[len(R)//2:]
        etacut = eta[len(eta)//2:]

        return np.mean(get_winding_numbers(Rcut, etacut, times, omega_values[omegaidx]))

    else:
        return 0


def get_phase_diagram(cc_directory, function_name, *args):

    """
    Obtain the values in the E0-omega phase diagram for a specific function.

    Args:
        cc_directory (str): The directory in which the collective coordinate data is contained
        function_name (function): The function for which the phase diagram should be calculated (e.g. get_mean_energy, max_radius)

    Returns:
        The 2D array of the phase diagram to be plotted

    """

    E0_values, omega_values = get_E0_and_omega_arrays(cc_directory + '/R/')

    phase_diagram_array = np.zeros((len(E0_values), len(omega_values)), dtype=float)

    for Eidx in range(len(E0_values)):
        print(function_name.__name__ + ' progress:', Eidx/len(E0_values), end='\r')
        for omegaidx in range(len(omega_values)):
            phase_diagram_array[Eidx, omegaidx] = function_name(cc_directory, E0_values, omega_values, Eidx, omegaidx, *args)

    return phase_diagram_array


def get_helicity_criterion(cc_directory, mean_winding_number_array):

    """
    Array which shows whether helicity is increasing on average, decreasing on average,
    and when it is monotonically or only on average increasing or decreasing.

    Args:
        cc_directory (str): The directory in which the collective coordinate data is contained
        mean_winding_number_array (ndarray): The array of mean winding numbers, output by
            mean_winding_number_array = get_phase_diagram(data_directory, mean_winding_number)

    Returns:
        2D array of the helicity criterion, used for plotting the phase diagram

    """

    E0_values, omega_values = get_E0_and_omega_arrays(cc_directory + '/R/')

    criterionarray = np.ones((len(E0_values), len(omega_values), 4), dtype=float)

    cmap = matplotlib.cm.get_cmap('PiYG')
    norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)

    monotonic_anticlockwise_colour = cmap(norm(1.))
    nonmonotonic_anticlockwise_colour = cmap(norm(.5))
    nonmonotonic_clockwise_colour = cmap(norm(-.5))
    monotonic_clockwise_colour = cmap(norm(-1.))

    criterionarray[np.where((mean_winding_number_array > 0.9))] = nonmonotonic_anticlockwise_colour
    criterionarray[np.where((mean_winding_number_array < -0.9))] = nonmonotonic_clockwise_colour

    for Eidx in range(len(E0_values)):
        print('Helicity criterion progress:', Eidx/len(E0_values), end='\r')
        for omegaidx in range(len(omega_values)):
            
            R, eta = get_coordinates(E0_values[Eidx], omega_values[omegaidx], cc_directory)

            etacut = eta[len(eta)//2:]

            # Monotonic changes
            if np.all(np.diff(etacut) > 0) and mean_winding_number_array[Eidx, omegaidx] > 0.9:
                criterionarray[Eidx, omegaidx] = monotonic_anticlockwise_colour
            elif np.all(np.diff(etacut) < 0) and mean_winding_number_array[Eidx, omegaidx] < -0.9:
                criterionarray[Eidx, omegaidx] = monotonic_clockwise_colour

    return criterionarray


def read_parameters(cc_directory):

    """
    Return various parameters related to the collective coordinate evolution such as external magnetic field
    and fitting parameters, stored, for the collective coordinate data, General/CCParams.txt, which is created
    from CCIntegration.sage.

    Args:
        cc_directory (str): The directory in which the collective coordinate data is contained

    Returns:
        All the parameters contained in CCParams.txt

    """

    # Read parameters from file
    with open(cc_directory + '/General/CCParams.txt') as f:

        for line in f.readlines():

            if 'Bz=' in line:
                Bz = float(line.split('=')[1])

            elif 'dw_width=' in line:
                dw_width = float(line.split('=')[1])

            elif 'E_exchangeFitParam' in line:
                E_exchangeFitParams = [float(param) for param in line.split('=')[1].split(',')]

            elif 'E_magnetic_integralFitParams' in line:
                E_magnetic_integralFitParams = [float(param) for param in line.split('=')[1].split(',')]

            elif 'E_electric_integralFitParams' in line:
                E_electric_integralFitParams = [float(param) for param in line.split('=')[1].split(',')]

            elif 'alpha' in line:
                alpha = float(line.split('=')[1])

            elif 'F_RexFitParams' in line:
                F_RexFitParams = [float(param) for param in line.split('=')[1].split(',')]

            elif 'Gamma11FitParams' in line:
                Gamma11FitParams = [float(param) for param in line.split('=')[1].split(',')]

            elif 'Gamma22FitParams' in line:
                Gamma22FitParams = [float(param) for param in line.split('=')[1].split(',')]

            elif 'G12FitParams' in line:
                G12FitParams = [float(param) for param in line.split('=')[1].split(',')]

            elif 'Xi1FitParams' in line:
                Xi1FitParams = [float(param) for param in line.split('=')[1].split(',')]

            elif 'Xi2FitParams' in line:
                Xi2FitParams = [float(param) for param in line.split('=')[1].split(',')]

            elif 'Xi3FitParams' in line:
                Xi3FitParams = [float(param) for param in line.split('=')[1].split(',')]

            elif 'Xi4FitParams' in line:
                Xi4FitParams = [float(param) for param in line.split('=')[1].split(',')]

            elif 'OptimalSkyrmionRadius' in line:
                optimalSkymionRadius = [float(param) for param in line.split('=')[1].split(',')]
                
    return Bz, dw_width, E_exchangeFitParams, E_magnetic_integralFitParams, E_electric_integralFitParams, \
        alpha, F_RexFitParams, Gamma11FitParams, Gamma22FitParams, G12FitParams, Xi1FitParams, Xi2FitParams, Xi3FitParams, Xi4FitParams, \
        optimalSkymionRadius
