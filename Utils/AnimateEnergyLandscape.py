"""
Title: Spacetime magnetic hopfions: from internal excitations and braiding of skyrmions

Authors and Affiliations: Ross Knapman[1,2,3], Timon Tausendpfund[1], Sebastián A. Díaz[2,4], Karin Everschor-Sitte[2,3]
    1 Institute of Physics, Johannes Gutenberg University Mainz, 55128 Mainz, Germany
    2 Faculty of Physics, University of Duisburg-Essen, 47057 Duisburg, Germany
    3 Center for Nanointegration Duisburg-Essen (CENIDE), University of Duisburg-Essen, 47057 Duisburg, Germany
    4 Department of Physics, University of Konstanz, 78457 Konstanz, Germany

Purpose: Render movie of energy landscape over time using collective coordinate modelling. Designed to be called by other
scripts rather than run directly.

Written 2023 by Ross Knapman <ross.knapman@uni-due.de>
"""

import numpy as np
from matplotlib import pyplot as plt
from Utils.Fits import *
import os
from multiprocessing import Pool


class Energy_Landscape_Parallel_Animator:

    """
    Class to create the energy landscape animation in parallel using multiprocessing.
    """

    def __init__(self, R_sol, eta_sol, Ez_array, Bz_array, start_idx, end_idx,
E_exchangeFitParams, E_magnetic_integralFitParams, E_electric_integralFitParams, out_name='EnergyLandscape.mp4'):

        self.R_sol = R_sol
        self.eta_sol = eta_sol
        self.Ez_array = Ez_array
        self.Bz_array = Bz_array
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.E_exchangeFitParams = E_exchangeFitParams
        self.E_magnetic_integralFitParams = E_magnetic_integralFitParams
        self.E_electric_integralFitParams = E_electric_integralFitParams
        self.out_name = out_name

        px = 1/plt.rcParams['figure.dpi']  # pixel in inches
        self.fig = plt.figure(figsize=(400*px, 400*px))
        self.ax = self.fig.add_subplot(projection='3d', azim=-90, elev=40) 

        self.R_3d = np.linspace(0.3, 1.2*np.max(R_sol[self.start_idx:self.end_idx]), 100)
        self.eta_3d = np.linspace(0., 2*np.pi, 100)
        self.R, self.ETA = np.meshgrid(self.R_3d, self.eta_3d)

        # Express the mesh in the cartesian system
        self.X, self.Y = self.R*np.cos(self.ETA), self.R*np.sin(self.ETA)

        # Get the minimum and maximum energy values (for choosing z-axis limits)
        self.minE = np.min(self._get_total_energy(np.max(Ez_array), np.min(Bz_array)))
        self.maxE = np.max(self._get_total_energy(np.min(Ez_array), np.max(Bz_array)))

        # Clear any clutter from previous renders
        os.system('rm /tmp/*.png')

    def _get_total_energy(self, Ez, Bz):

        exchange_energy = inverse_quadratic_linear_fit(self.R, self.E_exchangeFitParams[0], self.E_exchangeFitParams[1], self.E_exchangeFitParams[2], self.E_exchangeFitParams[3])
        magnetic_energy = -Bz * quadratic_fit(self.R, self.E_magnetic_integralFitParams[0], self.E_magnetic_integralFitParams[1], self.E_magnetic_integralFitParams[2])
        electric_energy = -Ez * np.cos(self.ETA) * linear_fit(self.R, self.E_electric_integralFitParams[0], self.E_electric_integralFitParams[1])
        return exchange_energy + magnetic_energy + electric_energy

    def _render_frame(self, i):

        self.ax.clear()
        R_idx = np.argmin(np.abs(self.R_sol[i] - self.R_3d))
        eta_idx = np.argmin(np.abs(self.eta_sol[i] % (2*np.pi) - self.eta_3d))
        energy_3d = self._get_total_energy(self.Ez_array[i], self.Bz_array[i])
        energy = energy_3d[eta_idx, R_idx]
        self.ax.plot_surface(self.X, self.Y, energy_3d, cmap=plt.cm.YlOrBr_r, alpha=0.6)
        self.ax.plot3D(self.R_sol[i]*np.cos(self.eta_sol[i]), self.R_sol[i]*np.sin(self.eta_sol[i]), energy, 'o', color='lime')

        self.ax.set_zlim(self.minE, self.maxE)
        self.ax.set_xlabel(r'$R\cos\eta$')
        self.ax.set_ylabel(r'$R\sin\eta$')
        self.ax.set_zlabel(r'$U$')

        plt.tight_layout()

        # Save to a temporary file
        plt.savefig('/tmp/%04d.png' %(i-self.start_idx))

    def render(self, no_processes=10):

        """
        Render frames of movie (in parallel using a multiprocessing pool) as individual images, then call ffmpeg to turn them into a movie.
        """

        with Pool(no_processes) as p:
            p.map(self._render_frame, [i for i in range(self.start_idx, self.end_idx)])

        os.system('ffmpeg -y -framerate 25 -i /tmp/%04d.png -r 25 -pix_fmt yuv420p ' + self.out_name + ' &>/dev/null')
