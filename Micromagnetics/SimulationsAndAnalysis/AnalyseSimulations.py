"""
Title: Spacetime magnetic hopfions: from internal excitations and braiding of skyrmions

Authors and Affiliations: Ross Knapman[1,2,3], Timon Tausendpfund[1], Sebastián A. Díaz[2,4], Karin Everschor-Sitte[2,3]
    1 Institute of Physics, Johannes Gutenberg University Mainz, 55128 Mainz, Germany
    2 Faculty of Physics, University of Duisburg-Essen, 47057 Duisburg, Germany
    3 Center for Nanointegration Duisburg-Essen (CENIDE), University of Duisburg-Essen, 47057 Duisburg, Germany
    4 Department of Physics, University of Konstanz, 78457 Konstanz, Germany

Purpose: Analysis of the micromagnetic simulations, called by RunSimulations.sh after running the simulation.
Outputs array of radius, helicity, and a movie of the micromagnetic simulation in the 'AnalysisResults' directory
of the corresponding simulation.

Written 2023 by Ross Knapman <ross.knapman@uni-due.de>
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
import discretisedfield as df
from scipy.constants import physical_constants
import os
import sys

sys.path.append('../..')
from Utils.QuantityReader import read_quantity

sys.path.append('../../Micromagnetics/SimulationsAndAnalysis')
from GeneralMicromagneticAnalysisTools import Read, Animate, Calculate

# Take in directory for which to analyse the data
simulation_dir = sys.argv[1]
params_file = '../../MicromagneticParams'

# In the simulation directory, create a subdirectory AnalysisResults (if it does not already exist)
try:
    os.mkdir(simulation_dir + '/AnalysisResults')
except FileExistsError:
    pass


##################################
# Get Parameters from Simulation #
##################################

# Get quantities to convert to dimensionless units
Ms = read_quantity('Ms', params_file)
I1 = read_quantity('I1', params_file)
I2 = read_quantity('I2', params_file)
gamma = physical_constants[u'electron gyromag. ratio'][0]
time_multiple = Ms * I2 / (gamma*I1*I1)
field_multiple = I1**2 / (Ms * I2)
length_multiple = np.sqrt(I2 / I1)
energy_multiple = np.sqrt(I1 * I2)


#############
# Load Data #
#############

directory = simulation_dir + '/Data/'
dx, dy, dz = Read.sampleDiscretisation(directory)
files_to_scan = Read.getFilesToScan(directory)
times = Read.simulationTimeArray(directory, loadMethod='files') / time_multiple
np.save(simulation_dir + '/AnalysisResults/times', times)


###############################
# Calculate Radius Variations #
###############################

radii = np.zeros(len(times), dtype=float)

for i in range(len(files_to_scan)):

    mz = df.Field.fromfile(directory + files_to_scan[i]).array[:, :, 0, 2]

    # Find contours using z-component of magnetization being zero (there should only be one such contour for a single skyrmion)
    contour = np.array(measure.find_contours(mz, level=0.))[0]

    # Get centre of skyrmion from which to calculate radius
    centroid_x = np.median(contour[:, 0])
    centroid_y = np.median(contour[:, 1])

    # Get distance between origin and points on contour
    x_displacements = (contour[:, 0] - centroid_x) * dx
    y_displacements = (contour[:, 1] - centroid_y) * dy

    contour_radii = np.sqrt(x_displacements**2 + y_displacements**2)
    radii[i] = np.mean(contour_radii)

np.save(simulation_dir + '/AnalysisResults/Radii.npy', radii / length_multiple)


#################################
# Calculate Helicity Variations #
#################################

helicities = Calculate.skyrmion_helicity_array(simulation_dir + '/Data')
np.save(simulation_dir + '/AnalysisResults/Helicities', helicities)


######################################################
# Animation of Magnetization Texture from Simulation #
######################################################

mean_radius = np.mean(radii)

px = 1/plt.rcParams['figure.dpi']  # pixel in inches
fig, ax = plt.subplots(figsize=(400*px, 400*px))

animator = Animate.MagnetizationAnimator('magnetization', directory, fig=fig, ax=ax, length_units=length_multiple, interpolation='bicubic', 
                                    plot_quiver=True, step=2, quiver_colour=[1, 1, 1], quiver_scale=3, quiver_headwidth=10, quiver_headlength=20,
                                    quiver_headaxislength=20, start_file='m009000.ovf', end_file='m009999.ovf', out_name=simulation_dir + '/AnalysisResults/Animation.mp4')

# Scale/translate axis labels to centre at centre of skyrmion show length in dimensionless units
plot_extent = animator.colour_plot.get_extent()
centre_X = (plot_extent[1] - plot_extent[0]) / 2
centre_Y = (plot_extent[3] - plot_extent[2]) / 2

animator.colour_plot.set_extent([plot_extent[0] - centre_X, plot_extent[1] - centre_X, plot_extent[2] - centre_Y, plot_extent[3] - centre_Y])
animator.quiver_plot.XY[:, 0] = animator.quiver_plot.XY[:, 0] - centre_X
animator.quiver_plot.XY[:, 1] = animator.quiver_plot.XY[:, 1] - centre_Y

axis_ticks = [i for i in range(-4, 5, 2)]

ax.set_xticks(axis_ticks)
ax.set_yticks(axis_ticks)

ax.set_xlim(np.min(axis_ticks), np.max(axis_ticks))
ax.set_ylim(np.min(axis_ticks), np.max(axis_ticks))

animator.animate()
