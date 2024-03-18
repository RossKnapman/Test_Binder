"""
Title: Spacetime magnetic hopfions: from internal excitations and braiding of skyrmions

Authors and Affiliations: Ross Knapman[1,2,3], Timon Tausendpfund[1], Sebastián A. Díaz[2,4], Karin Everschor-Sitte[2,3]
    1 Institute of Physics, Johannes Gutenberg University Mainz, 55128 Mainz, Germany
    2 Faculty of Physics, University of Duisburg-Essen, 47057 Duisburg, Germany
    3 Center for Nanointegration Duisburg-Essen (CENIDE), University of Duisburg-Essen, 47057 Duisburg, Germany
    4 Department of Physics, University of Konstanz, 78457 Konstanz, Germany

Purpose: Renders movies of energy landscape and POV-Ray render of skyrmion motion from collective coordinate integration
data, which is put in a directory "Movies".

Written 2023 by Ross Knapman <ross.knapman@uni-due.de>
"""


import numpy as np
import multiprocessing
import os

import sys
sys.path.append('..')
from Utils.Fits import *
from Utils.AnimateEnergyLandscape import Energy_Landscape_Parallel_Animator
from Utils.PhaseDiagramFunctions import *
from Utils.POVRender import POVRenderer


if __name__ == '__main__':

    # The points (E0, omega) to be analysed
    parameter_pairs = [[0.25, 0.72], [0.20, 1.18]]

    # Number of processes used for rendering POV-Ray and energy surface videos
    no_processes = multiprocessing.cpu_count()

    # Directory in which videos for each point should be stored
    out_directory = 'Movies'
    os.makedirs(out_directory, exist_ok=True)

    # Directory in which data is contained
    data_directory = 'CollectiveCoordinateIntegration/CollectiveCoordinateDataSupplied'

    # Times array
    times = np.load(data_directory + '/General/times.npy')

    # Loop through the phase diagram points defined in parameter_pairs list above
    for i in range(len(parameter_pairs)):

        E0 = parameter_pairs[i][0]
        omega = parameter_pairs[i][1]

        # Obtain parameters such as magnetic field, fit parameters (see Methods section of manuscript on fitting)
        # E_exchange, E_magnetic, and E_electric are the contributions from exchange, Zeeman, and electric field to energy respectively
        # alpha is the Gilbert damping constant
        # F_rex is -dU_ex / dR
        # Gamma11 and Gamma22 are as in the manuscript where 1 is R and 2 is η
        # G12 is the G_{Rη} from the manuscript

        # What I call Xi are parameters that show up in the analytical expressions for the helicity-dependent
        # parts of the generalised forces for the skyrmion profile ansatz we use

        # Xi1 is the integral of cos(2θ) dθ / dR over ρ
        # Xi2 is the integral of ρ * d^2θ / dRdρ over ρ
        # Xi3 is the integral of cos(θ)sin(θ) over ρ
        # Xi4 is the integral of ρ * dθ / dρ over ρ
        Bz, dw_width, E_exchangeFitParams, E_magnetic_integralFitParams, E_electric_integralFitParams, alpha, \
            F_RexFitParams, Gamma11FitParams, Gamma22FitParams, G12FitParams, Xi1FitParams, Xi2FitParams, \
            Xi3FitParams, Xi4FitParams, optimalSkyrmionRadius = read_parameters(data_directory)

        # Create the directory for this point on the phase diagram
        out_directory_point = out_directory + '/E0%.2fOmega%.2f' %(E0, omega)
        if not os.path.exists(out_directory_point):
            os.mkdir(out_directory_point)

        # Arrays of magnetic and electric fields over collective coordinate integration
        Bz_array = np.full_like(times, Bz)
        Ez_array = E0 * np.cos(omega*times)

        # Load in collective coordinate arrays
        R, eta = get_coordinates(E0, omega, data_directory)

        # Calculate array of skyrmion energy during the collective coordinate integration
        energies = total_energy(R, eta, Bz_array, Ez_array, E_exchangeFitParams, E_magnetic_integralFitParams, E_electric_integralFitParams)


        ##############################################################
        # Animations of Energy Landscape and Skyrmion Motion Renders #
        ##############################################################

        # Animate for last 10% of times
        start_idx = int(0.9*len(R))
        end_idx = len(R)

        # Animate energy landscape
        energy_animator = Energy_Landscape_Parallel_Animator(R, eta, Ez_array, Bz_array, start_idx, end_idx, E_exchangeFitParams, E_magnetic_integralFitParams, E_electric_integralFitParams, out_name=out_directory_point+'/EnergyLandscape.mp4')
        energy_animator.render(no_processes=no_processes)

        # POV-Ray render of skyrmion motion
        povrenderer = POVRenderer(R, eta, times, dw_width, start_idx=start_idx, end_idx=end_idx, out_name=out_directory_point+'/POV.mp4')
        povrenderer.render()
