"""
Title: Spacetime magnetic hopfions: from internal excitations and braiding of skyrmions

Authors and Affiliations: Ross Knapman[1,2,3], Timon Tausendpfund[1], Sebastián A. Díaz[2,4], Karin Everschor-Sitte[2,3]
    1 Institute of Physics, Johannes Gutenberg University Mainz, 55128 Mainz, Germany
    2 Faculty of Physics, University of Duisburg-Essen, 47057 Duisburg, Germany
    3 Center for Nanointegration Duisburg-Essen (CENIDE), University of Duisburg-Essen, 47057 Duisburg, Germany
    4 Department of Physics, University of Konstanz, 78457 Konstanz, Germany

Purpose: Tools to create movie of skyrmion using collective coordinate data, rendered using POV-Ray.

Written 2023 by Ross Knapman <ross.knapman@uni-due.de>
"""

import numpy as np
import os
from multiprocessing import Pool
from colorsys import hls_to_rgb

class POVRenderer:

    """
    Class to create a movie of the skyrmion (with magnetization represented as cones) using the collective
    coordinate time evolution data using POV-Ray.
    """

    def __init__(self, R_sol, eta_sol, times, dw_width, no_rings=5, N1=5, outer_theta=0.1, start_idx=0, end_idx=None, out_name='POV.mp4'):
        self.R_sol = R_sol
        self.eta_sol = eta_sol
        self.times = times
        self.dw_width = dw_width
        self.no_rings = no_rings  # Number of rings of spins
        self.N1 = N1  # Number of vectors to put on a ring of radius 1 (i.e. sets how close together)
        self.outer_theta = outer_theta  # Theta which defines position of outer radius

        self.start_idx = start_idx
        self.end_idx = len(self.R_sol) if not end_idx else end_idx
        self.out_name = out_name

        # Calculate how many vectors on each radius ring, using the mean value of the radius
        self.mean_R = np.mean(self.R_sol)
        self.mean_outer_radius = self.dw_width * np.arcsinh(np.sinh(self.mean_R/self.dw_width) / np.tan(outer_theta / 2))
        self.mean_rho_values = np.arange(0, self.mean_outer_radius, self.mean_outer_radius / no_rings)
        self.N_rho = (self.mean_rho_values*N1).astype(int)

        # Calculate the angle required for the skymion to remain in the field of view, if camera at a height of 15
        self.all_outer_radii = self.dw_width * np.arcsinh(np.sinh(self.R_sol/self.dw_width) / np.tan(self.outer_theta / 2))
        self.view_angle = 2 * np.arctan(np.max(self.all_outer_radii[start_idx:end_idx]) / 15) * 180 / np.pi

        # Clear any clutter from previous renders
        os.system('rm /tmp/*.png')
        os.system('rm /tmp/*.pov')


    def render_frame(self, i):

        """
        Render the ith frame.
        """

        R = self.R_sol[i]
        eta = self.eta_sol[i]

        outer_radius = self.all_outer_radii[i]
        rho_values = np.arange(0, outer_radius, outer_radius / self.no_rings)

        with open('/tmp/%04d.pov' %i, 'w') as outfile:

            module_path = os.path.abspath(__file__)

            with open(os.path.dirname(module_path) + '/Template.pov') as infile:

                for line in infile.readlines():
                    if '{{ Angle }}' in line:
                        line_to_write = line.replace("{{ Angle }}", str(self.view_angle))
                        outfile.write(line_to_write)
                    
                    else:
                        outfile.write(line)

                for rho_idx in range(len(rho_values)):
                    
                    rho = rho_values[rho_idx]

                    # Get number of vectors at this radius
                    Nr = self.N_rho[rho_idx]

                    if Nr != 0:
                        angle_increment = 2*np.pi / Nr
                        psi_values = np.arange(0, 2*np.pi, angle_increment)
                    else:  # Centre spin
                        psi_values = np.array([0.])

                    for psi in psi_values:
                        x = rho*np.cos(psi)
                        y = rho*np.sin(psi)

                        if not np.isclose(rho, 0):
                            theta = 2*np.arctan(np.sinh(R/self.dw_width) / np.sinh(rho/self.dw_width))
                        else:
                            theta = np.pi
                        phi = psi + eta

                        # Get HLS values
                        s = 1.
                        l = 0.5*np.cos(theta) + 0.5
                        h = phi / (2*np.pi)

                        # Convert to RGB
                        r, g, b = hls_to_rgb(h, l, s)

                        # Append vector to render file
                        outfile.write('Vector(' \
                            + str(x) + ',' \
                            + str(y) + ',' \
                            + str(phi) + ',' \
                            + str(theta) + ',' \
                            + str(r) + ',' \
                            + str(g) + ',' \
                            + str(b) + ')\n')

        os.system('povray +O"/tmp/%04d.png" +UA +A +H400 +W400 /tmp/%04d.pov &>/dev/null' %(i, i))

        # Add a white background to the frame (requires ImageMagick)
        os.system('convert /tmp/%04d.png -background white -alpha remove /tmp/%04dWhiteBG.png &>/dev/null' %(i, i-self.start_idx))
        

    def _render_all_frames(self, parallel, no_processes):

        # Process in parallel
        if parallel:
            
            with Pool(no_processes) as p:
                p.map(self.render_frame, [i for i in range(self.start_idx, self.end_idx)])

        else:
            # Loop through frames, rendering each
            for i in range(self.start_idx, self.end_idx):
                print(i, end='\r')
                self.render_frame(i)


    def render(self, parallel=True, no_processes=10):

        """
        Render the movie after having created the individual frames (which are stored in /tmp/).
        """

        self._render_all_frames(parallel, no_processes)
        os.system('ffmpeg -y -framerate 25 -i /tmp/%04dWhiteBG.png -r 25 -pix_fmt yuv420p ' + self.out_name + ' &>/dev/null')
                    