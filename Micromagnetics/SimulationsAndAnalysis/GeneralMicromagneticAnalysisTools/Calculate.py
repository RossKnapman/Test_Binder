"""
Title: Spacetime magnetic hopfions: from internal excitations and braiding of skyrmions

Authors and Affiliations: Ross Knapman[1,2,3], Timon Tausendpfund[1], Sebastián A. Díaz[2,4], Karin Everschor-Sitte[2,3]
    1 Institute of Physics, Johannes Gutenberg University Mainz, 55128 Mainz, Germany
    2 Faculty of Physics, University of Duisburg-Essen, 47057 Duisburg, Germany
    3 Center for Nanointegration Duisburg-Essen (CENIDE), University of Duisburg-Essen, 47057 Duisburg, Germany
    4 Department of Physics, University of Konstanz, 78457 Konstanz, Germany

Purpose: Collection of functions to perform calculations on micromagnetic simulation data (such as skyrmion position).

Written 2023 by Ross Knapman <ross.knapman@uni-due.de>
"""

import numpy as np
import discretisedfield as df
import findiff

from GeneralMicromagneticAnalysisTools import Read


def skyrmionCOM(directory, inFile, dx, dy, zIndex=0, edgeCutXFrac=0.1, edgeCutYFrac=0.1):
    """Calculate the centre of the skyrmion by obtaining the "centre of mass" of the skyrmion number density.

    Args:
        directory (str): The directory in which the simulation data is stored.
        inFile (str): The filename for which to calculate the skyrmion's centre.
        dx (float64): The simulation cell size in the x-dimension, in nm.
        dy (float64): The simulation cell size in the y-dimension, in nm.
        zIndex (int, optional): The index along the z-axis for which the skyrmion centre should be calculated (for a 3D sample).
        edgeCutXFrac (float64, optional): The fraction of the system along the x-axis that should be cut off at the edges (thus allowing edge effects to be excluded).
        edgeCutYFrac (float64, optional): Same as above but for the y-axis.

    Returns:
        A two-element array of the form `[x-coordinate of centre of skyrmion, y-coordinate of centre of skyrmion]`.

    """

    Lx, Ly = Read.sampleExtent(directory)

    edgeCutX = int(np.round(edgeCutXFrac * Lx))
    edgeCutY = int(np.round(edgeCutYFrac * Ly))

    m = df.Field.fromfile(directory + '/' + inFile).array[:, :, zIndex]
    m = m[edgeCutX:m.shape[0]-edgeCutX, edgeCutY:m.shape[1]-edgeCutY]

    # Define operators
    # args: axis, lattice constant, derivative order, accuracy
    d_dx = findiff.FinDiff(0, 1, 1, acc=4)
    # args: axis, lattice constant, derivative order, accuracy
    d_dy = findiff.FinDiff(1, 1, 1, acc=4)

    # Apply them to m
    mdx = d_dx(m)
    mdy = d_dy(m)

    rhoSk = np.einsum('ijk,ijk->ij', m, np.cross(mdx, mdy))
    Nsk = np.sum(rhoSk)  # Out from the "usual" Nsk by 4pi

    xy = np.indices((m.shape[0], m.shape[1])).transpose(1, 2, 0)
    comArray = np.einsum('ijk,ij->ijk', xy, rhoSk)

    return np.array([(np.sum(comArray[:, :, 0] / Nsk) + edgeCutX) * dx * 1e9, (np.sum(comArray[:, :, 1] / Nsk) + edgeCutY) * dy * 1e9])


def skyrmion_helicity(directory, filename):

    """Calculate the skyrmion helicity from an ovf file.
    The calculation works by sweeping down the y-axis, and getting the points to the left and right of the
    centre, which are closest to m_z = 0 (where we define the radius of the skyrmion). We then take the mean
    helicity of all of these points on the radius. We assume that the system is a thin film, i.e. that the
    system only has one layer along the z-axis.

    IMPORTANT: We assume a single skyrmion in a collinear background!
    
    Args:
        directory (str): The directory containing the ovf file for which the helicity should be calculated.
        filename (str): The ovf file containing the skyrmion for which the helicity should be calculated.

    Returns:
        The calculated helicity (between -pi and pi).

    """

    # Load the magnetization array
    m = df.Field.fromfile(directory + '/' + filename).array

    # Read the cell size
    dx, dy, dz = Read.sampleDiscretisation(directory)

    # Calculate the centre of mass of the skyrmion
    com = skyrmionCOM(directory, filename, dx, dy, edgeCutXFrac=0., edgeCutYFrac=0.)

    # Get the indices of the cells closest to the centre of mass
    central_x_idx = int(np.round(com[0] / (dx*1e9)))
    central_y_idx = int(np.round(com[1] / (dy*1e9)))

    # Get the core polarization of the skyrmion
    polarization = -1 if m[central_x_idx, central_y_idx, 0, 2] < 0 else 1

    # List to store points on the radius (as indices of array)
    radius_points = []

    # Loop through y-values, getting x-values that are on the skyrmion boundary (thanks to Robin Msiska for the inspration)
    for y_idx in range(m.shape[1]):
    
        # Check that there are parts of this line that are actually within the skyrmion radius (defined by m_z = 0)
        if np.any(np.sign(m[:, y_idx, 0, 2]) == np.sign(polarization)):

            # Get points to left and right of centre of skyrmion, on its radius
            radius_points.append([np.argmin(np.abs(m[:central_x_idx, y_idx, 0, 2])), y_idx])
            radius_points.append([central_x_idx + np.argmin(np.abs(m[central_x_idx:, y_idx, 0, 2])), y_idx])

    # Get helicities for each point
    helicities = np.zeros(len(radius_points), dtype=float)

    for i in range(len(helicities)):

        x_idx = radius_points[i][0]
        y_idx = radius_points[i][1]

        # Get x- and y-position of point in nm
        x = x_idx * dx * 1e9
        y = y_idx * dy * 1e9

        # Get displacements from centre of skyrmion
        x_displacement = x - com[0]
        y_displacement = y - com[1]

        # Calculate polar coordinate in plane
        plane_angle = np.arctan2(y_displacement, x_displacement)
        
        # Calculate polar coordinate of in-plane components of magnetization
        phi = np.arctan2(m[x_idx, y_idx, 0, 1], m[x_idx, y_idx, 0, 0])

        # Save helicity values
        helicities[i] = phi - plane_angle

    # Deal with helicities that are outside of the range [-pi, pi]
    for i in range(len(helicities)):
        while helicities[i] > np.pi:
            helicities[i] -= 2*np.pi            
        while helicities[i] < -np.pi:
            helicities[i] += 2*np.pi

    # Below, deal with the case that the helicity is close to pi, such that some values have
    # helicity ≈ -pi and other have helicity ≈ pi

    # If absolute value of helicity is close to pi, and there is a mixture of positive and negative helicities
    if np.abs(np.mean(np.abs(helicities)) - np.pi) < np.pi/4 and np.any(helicities > 0) and np.any(helicities < 0):

        # Instead get angle w.r.t. pi direction (-x), to avoid discontinous point at helicity = pi or -pi
        anglesToMinusX = np.zeros_like(helicities, dtype=float)

        for i in range(len(anglesToMinusX)):

            if helicities[i] >= 0:
                anglesToMinusX[i] = np.pi - helicities[i]
            
            else:
                anglesToMinusX[i] = -np.pi - helicities[i]
                
        averageAngleToMinusX = np.mean(anglesToMinusX)

        if averageAngleToMinusX >= 0:
            return np.pi - averageAngleToMinusX
        
        else:
            return -np.pi - averageAngleToMinusX

    else:
        return np.average(helicities)


def skyrmion_helicity_array(directory, startFile=None, endFile=None):

    """Get an array of skyrmion helicities for a given simulation directory.
    Note that this only works for a single skyrmion in a collinear background.

    Args:
        directory (str): The directory containing the ovf file for which the helicity should be calculated.
        startFile (str, optional): The starting file for which the helicity should be calculated.
        endFile (str, optional): The ending file for which the helicity should be calculated.

    Returns:
        Array of calculated helicities.

    """

    filesToScan = Read.getFilesToScan(directory, startFile, endFile)

    helicities = np.zeros(len(filesToScan), dtype=float)

    for i in range(len(filesToScan)):

        print("Calculating helicity", i, "of", len(filesToScan) - 1, end='\r')
        helicities[i] = skyrmion_helicity(directory, filesToScan[i])

    return helicities
