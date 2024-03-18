"""
Title: Spacetime magnetic hopfions: from internal excitations and braiding of skyrmions

Authors and Affiliations: Ross Knapman[1,2,3], Timon Tausendpfund[1], Sebastián A. Díaz[2,4], Karin Everschor-Sitte[2,3]
    1 Institute of Physics, Johannes Gutenberg University Mainz, 55128 Mainz, Germany
    2 Faculty of Physics, University of Duisburg-Essen, 47057 Duisburg, Germany
    3 Center for Nanointegration Duisburg-Essen (CENIDE), University of Duisburg-Essen, 47057 Duisburg, Germany
    4 Department of Physics, University of Konstanz, 78457 Konstanz, Germany

Purpose: Collection of functions used in the Jupyter notebooks that create the manuscript figures.

Written 2023 by Ross Knapman <ross.knapman@uni-due.de>
"""

import numpy as np
from Utils.PhaseDiagramFunctions import *


def plot_evolution(cc_directory, mumax_directory, ax_evolution, E0, omega, min_aspect_ratio, plot_colour):

    """
    Plot Rcos(eta)-Rsin(eta).

    Args:
        cc_directory (str): The directory in which the collective coordinate data is contained
        mumax_directory (str): The directory in which the micromagnetic data is contained
        ax_evolution (str): The axis on which the plot of Rcos(eta)-Rsin(eta) should be plotted
        E0 (float): The electric field amplitude
        omega (float): The electric field angular frequency
        min_aspect_ratio (float): The minimum value of max(abs(Rcos(eta))) / max(abs(Rsin(eta))) for all the points to be plotted
        plot_colour (str): The colour of line of the plot

    """

    R_cos_eta_cc, R_sin_eta_cc, R_cos_eta_mumax, R_sin_eta_mumax = get_cut_coordinate_arrays(cc_directory, mumax_directory, E0, omega)

    ax_evolution.plot(R_cos_eta_cc, R_sin_eta_cc, '--', color=plot_colour, label='Collective Coordinates')
    ax_evolution.plot(R_cos_eta_mumax, R_sin_eta_mumax, '-', color=plot_colour, linewidth=3, label='Micromagnetics')

    max_xvalue, max_yvalue = get_evolution_plot_axis_limits(cc_directory, mumax_directory, E0, omega)

    # As we want to have equal aspect ratio in plots, but also want the plots to align, we adjust the
    # y-axis limit so that the plots have the same height, whilst preserving the aspect ratio, so we take
    # the minimum aspect ratio from all the plots
    aspect_ratio = (max_xvalue / max_yvalue) / min_aspect_ratio

    # Only integer values of xticks (I just do -5 to 5 as I know this range is bigger than what will actually be plotted)
    ax_evolution.set_xticks([i for i in range(-5, 6)])
    ax_evolution.set_yticks([i for i in range(-5, 6)])

    ax_evolution.set_xlim(-1.1*max_xvalue, 1.1*max_xvalue)
    ax_evolution.set_ylim(-1.1*aspect_ratio * max_yvalue, 1.1*aspect_ratio * max_yvalue)

    ax_evolution.set_xlabel(r'$R\cos\eta$')
    ax_evolution.set_ylabel(r'$R\sin\eta$')
    
    ax_evolution.set_aspect('equal')


def get_start_idx(eta_sol, times, omega):

    """
    Take where helicity is closest to zero, but time is greater than 500 and less than 1000 - (time period of E_z).

    Args:
        eta_sol (ndarray): Array of skyrmion helicities obtained by collective coordinate integration
        times (ndarray): Array of time of the simulation
        omega (float): Angular frequency of the electric field

    Returns:
        The index of the simulation starting point within the (potentially cropped) array of helicity solutions.

    """
    return np.argmin(np.abs(eta_sol[np.where((times > times[len(times)//2]) & (times < np.max(times) - 2*np.pi/omega))] % (2*np.pi)))


def get_coordinate_plot_data(cc_directory, mumax_directory, E0, omega):

    """
    Loads the collective coordinate data for a given E0, omega.

    Args:
        cc_directory (str): The directory in which the collective coordinate data is contained
        mumax_directory (str): The directory in which the micromagnetic data is contained
        E0 (float): The electric field amplitude
        omega (float): The electric field angular frequency

    Returns:
        The array of skyrmion R from micromagnetics, skyrmion helicity from micromagnetics,
        skyrmion R from collective coordinates, skyrmion helicity from collective coordinates.

    """

    directory = mumax_directory + '/E0%.2fomega%.2f' %(E0, omega)
    radii_mumax = np.load(directory + '/AnalysisResults/Radii.npy')
    helicities_mumax = np.load(directory + '/AnalysisResults/Helicities.npy')
    radii_cc, helicities_cc = get_coordinates(float(E0), float(omega), cc_directory)

    return radii_mumax, helicities_mumax, radii_cc, helicities_cc


def get_cut_coordinate_arrays(cc_directory, mumax_directory, E0, omega):

    """
    Obtains R and eta from both CC and micromagnetic data and cuts it down to one period of the electric field.

    Args:
        cc_directory (str): The directory in which the collective coordinate data is contained
        mumax_directory (str): The directory in which the micromagnetic data is contained
        E0 (float): The electric field amplitude
        omega (float): The electric field angular frequency

    Returns:
        Array of R*cos(eta) from collective coordinates, R*sin(eta) from collecive coordinates,
        R*cos(eta) from micromagnetics, R*sin(eta) from micromagnetics

    """

    radii_mumax, helicities_mumax, radii_cc, helicities_cc = get_coordinate_plot_data(cc_directory, mumax_directory, E0, omega)

    times = np.load(cc_directory + '/General/times.npy')

    indices_per_time_interval = int(np.round(1/(times[1]-times[0])))
    start_idx = get_start_idx(helicities_cc, times, omega)
    time_start = times[start_idx] + times[len(times)//2]

    # Update the start index to the "correct value", after 1/2 of the integration time
    start_idx = int(np.round(time_start * indices_per_time_interval))

    # End index is just start index plus one period of the electric field oscillation
    end_idx = int(np.round((time_start + 2*np.pi/omega)*indices_per_time_interval)) + 1

    # Trim the quantity arrays (from both collective coordinates and micromagnetics) to between the start and end indices
    radii_cc = radii_cc[start_idx:end_idx]
    helicities_cc = helicities_cc[start_idx:end_idx]
    radii_mumax = radii_mumax[start_idx:end_idx]
    helicities_mumax = helicities_mumax[start_idx:end_idx]

    R_cos_eta_cc = radii_cc * np.cos(helicities_cc)
    R_sin_eta_cc = radii_cc * np.sin(helicities_cc)
    R_cos_eta_mumax = radii_mumax * np.cos(helicities_mumax)
    R_sin_eta_mumax = radii_mumax * np.sin(helicities_mumax)

    return R_cos_eta_cc, R_sin_eta_cc, R_cos_eta_mumax, R_sin_eta_mumax


def get_evolution_plot_axis_limits(cc_directory, mumax_directory, E0, omega):

    """
    Obtains the minimum and maximum values in x and y for the evolution plots.
    Return Rcos(eta) and Rsin(eta) for both CC and micromagnetics.

    Args:
        cc_directory (str): The directory in which the collective coordinate data is contained
        mumax_directory (str): The directory in which the micromagnetic data is contained
        E0 (float): The electric field amplitude
        omega (float): The electric field angular frequency

    Returns:
    The minimum and maximum (absolute) values of R*cos(eta) and R*sin(eta) considering both micromagnetics
    and collective coordinates

    """

    R_cos_eta_cc, R_sin_eta_cc, R_cos_eta_mumax, R_sin_eta_mumax = get_cut_coordinate_arrays(cc_directory, mumax_directory, E0, omega)

    max_xvalue = np.maximum(np.max(np.abs(R_cos_eta_cc)), np.max(np.abs(R_cos_eta_mumax)))
    max_yvalue = np.maximum(np.max(np.abs(R_sin_eta_cc)), np.max(np.abs(R_sin_eta_mumax)))

    return max_xvalue, max_yvalue

