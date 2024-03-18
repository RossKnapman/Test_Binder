"""
Title: Spacetime magnetic hopfions: from internal excitations and braiding of skyrmions

Authors and Affiliations: Ross Knapman[1,2,3], Timon Tausendpfund[1], Sebastián A. Díaz[2,4], Karin Everschor-Sitte[2,3]
    1 Institute of Physics, Johannes Gutenberg University Mainz, 55128 Mainz, Germany
    2 Faculty of Physics, University of Duisburg-Essen, 47057 Duisburg, Germany
    3 Center for Nanointegration Duisburg-Essen (CENIDE), University of Duisburg-Essen, 47057 Duisburg, Germany
    4 Department of Physics, University of Konstanz, 78457 Konstanz, Germany

Purpose: SageMath script to run the collective coordinate integrations in parallel. Calls a function defined in Fit.c, which must
be compiled with
    gcc -fPIC -shared -o c_fit.o Fit.c

Written 2023 by Ross Knapman <ross.knapman@uni-due.de>
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp
from ctypes import c_double, cdll
from numpy.ctypeslib import ndpointer
import multiprocessing
from multiprocessing import Pool
import os

import sys
sys.path.append('../../')
from Utils.Fits import *


# We use a C implementation of the time_derivatives function (called by ODE solver) to speed up the integration
# Compile using gcc -fPIC -shared -o c_fit.o Fit.c
lib = cdll.LoadLibrary('./c_fit.o')
CTimeDerivatives = lib.time_derivatives
CTimeDerivatives.restype = ndpointer(dtype=c_double, shape=(2,))  # Result type


def get_thiele(Bz, E0, alpha, omega):

    #########################################################
    # Use Total Skyrmion Energy to Calculate Optimal Radius #
    #########################################################

    # Finer mesh of value of skyrmion radii to more accurately obtain that which minimises energy
    R_finer = np.linspace(np.min(R_values), np.max(R_values), 1000)

    # Arrays of various contributions to energies (obtained from the fitting run at the start)
    E_exchange_values = inverse_quadratic_linear_fit(R_finer, E_exchangeFitParams[0], E_exchangeFitParams[1], E_exchangeFitParams[2], E_exchangeFitParams[3])
    E_magnetic_values = -Bz * quadratic_fit(R_finer, E_magnetic_integralFitParams[0], E_magnetic_integralFitParams[1], E_magnetic_integralFitParams[2])
    E_electric_values = -E0 * linear_fit(R_finer, E_electric_integralFitParams[0], E_electric_integralFitParams[1])
    E_values = E_exchange_values + E_magnetic_values - E_electric_values
    
    
    def time_derivatives(t, y):
    
        """
        Wrapper for C function call; calculates dR_dt and dη_dt.
        """

        R, eta = y

        Bz = Bz_array[np.argmin(np.abs(t - times))]
        Ez = Ez_array[np.argmin(np.abs(t - times))]

        return CTimeDerivatives(c_double(t), c_double(R), c_double(eta), c_double(alpha), c_double(Ez), c_double(Bz),
                               c_double(Gamma11FitParams[0]), c_double(Gamma11FitParams[1]), c_double(Gamma11FitParams[2]), c_double(Gamma11FitParams[3]),
                               c_double(Gamma22FitParams[0]), c_double(Gamma22FitParams[1]),
                               c_double(G12FitParams[0]), c_double(G12FitParams[1]),
                               c_double(F_RexFitParams[0]), c_double(F_RexFitParams[1]), c_double(F_RexFitParams[2]), c_double(F_RexFitParams[3]), c_double(F_RexFitParams[4]), c_double(F_RexFitParams[5]),
                               c_double(Xi1FitParams[0]), c_double(Xi1FitParams[1]), c_double(Xi1FitParams[2]), c_double(Xi1FitParams[3]), c_double(Xi1FitParams[4]), c_double(Xi1FitParams[5]),
                               c_double(Xi2FitParams[0]), c_double(Xi2FitParams[1]), c_double(Xi2FitParams[2]), c_double(Xi2FitParams[3]), c_double(Xi2FitParams[4]), c_double(Xi2FitParams[5]),
                               c_double(Xi3FitParams[0]), c_double(Xi3FitParams[1]), c_double(Xi3FitParams[2]), c_double(Xi3FitParams[3]), c_double(Xi3FitParams[4]), c_double(Xi3FitParams[5]),
                               c_double(Xi4FitParams[0]), c_double(Xi4FitParams[1]))

    optimal_R = R_finer[np.argmin(E_values)]  # Radius that minimises energy
    init = [optimal_R, np.pi]                 # Initial R and eta

    np.save(over_directory + '/General/times', times)

    # Magnetic field array (constant Bz in this case)
    Bz_array = np.zeros_like(times, dtype=float)
    Bz_array[:] = Bz

    # Electric field array
    Ez_array = np.zeros_like(times, dtype=float)
    Ez_array[:] = E0 * cos(omega*times)

    # Perform the time integration
    sol = solve_ivp(time_derivatives, [np.min(times), np.max(times)], init, t_eval=times, method='Radau')
    
    return sol


def initialiser(over_directory_in, times_in, omega_values_in, Bz_in, E0_values_in, alpha_in, R_values_in,
    E_exchangeFitParams_in, E_magnetic_integralFitParams_in, E_electric_integralFitParams_in,
    Gamma11FitParams_in, Gamma22FitParams_in, G12FitParams_in, F_RexFitParams_in,
    Xi1FitParams_in, Xi2FitParams_in, Xi3FitParams_in, Xi4FitParams_in):

    """
    Called by each multiprocessing.Pool worker when it starts (to initialise times, fields, fit parameters etc. for each thread).
    See comment under section "Fit Quantities" for explanation of the various values.
    """

    global over_directory
    global times
    global omega_values
    global Bz
    global E0_values
    global alpha
    global R_values
    global E_exchangeFitParams
    global E_magnetic_integralFitParams
    global E_electric_integralFitParams
    global Gamma11FitParams
    global Gamma22FitParams
    global G12FitParams
    global F_RexFitParams
    global Xi1FitParams
    global Xi2FitParams
    global Xi3FitParams
    global Xi4FitParams

    over_directory = over_directory_in
    times = times_in
    omega_values = omega_values_in
    Bz = Bz_in
    E0_values = E0_values_in
    alpha = alpha_in
    R_values = R_values_in
    E_exchangeFitParams = E_exchangeFitParams_in
    E_magnetic_integralFitParams = E_magnetic_integralFitParams_in
    E_electric_integralFitParams = E_electric_integralFitParams_in
    Gamma11FitParams = Gamma11FitParams_in
    Gamma22FitParams = Gamma22FitParams_in
    G12FitParams = G12FitParams_in
    F_RexFitParams = F_RexFitParams_in
    Xi1FitParams = Xi1FitParams_in
    Xi2FitParams = Xi2FitParams_in
    Xi3FitParams = Xi3FitParams_in
    Xi4FitParams = Xi4FitParams_in


def evolve(Eidx):

    """
    For the electric field input, loop through all omega values and perform the integration.
    This is the function called by each multiprocessing Pool process.
    """

    for omegaidx in range(len(omega_values)):

        filename = 'E0{:.2f}'.format(E0_values[Eidx]) + 'omega{:.2f}'.format(omega_values[omegaidx])

        if os.path.exists(over_directory + '/R/' + filename + '.npy') and os.path.exists(over_directory + '/eta/' + filename + '.npy'):
            print('Already exists', filename)

        else:
            print('Does not exist', filename)
            sol = get_thiele(Bz, E0_values[Eidx], alpha, omega_values[omegaidx])
            np.save(over_directory + '/R/E0{:.2f}'.format(E0_values[Eidx]) + 'omega{:.2f}'.format(omega_values[omegaidx]), sol.y[0, :])
            np.save(over_directory + '/eta/E0{:.2f}'.format(E0_values[Eidx]) + 'omega{:.2f}'.format(omega_values[omegaidx]), sol.y[1, :])


if __name__ == '__main__':

    over_directory = 'CollectiveCoordinateData'

    if not os.path.exists(over_directory):

        os.mkdir(over_directory)
        os.mkdir(over_directory + '/R')
        os.mkdir(over_directory + '/eta')
        os.mkdir(over_directory + '/General')


    ##################################################################
    # Define Variables for Symbolic Algebra Required for Integration #
    ##################################################################

    # Here, the variables are defined, such as polar coordinate variables, helicity etc.

    print('Defining variables')

    rho = var('rho')  # Polar coordinate radius
    psi = var('psi')  # Polar coordinate angle

    m   = var('m')    # Skyrmion vorticity
    eta = var('eta')  # Skyrmion helicity
    R   = var('R')    # Skyrmion radius
    w   = var('w')    # Skyrmion domain wall width

    # Theta and phi as in manuscript
    theta = 2 * arctan(sinh(R/w) / sinh(rho/w))
    phi = m*psi + eta

    # Components of magnetization
    mx = cos(phi) * sin(theta)
    my = sin(phi) * sin(theta)
    mz = cos(theta)

    # Components of Laplacian in polar coordinates
    mxLaplacian = diff(mx, rho, 2) + (1/rho) * diff(mx, rho) + (1/rho^2) * diff(mx, psi, 2)
    myLaplacian = diff(my, rho, 2) + (1/rho) * diff(my, rho) + (1/rho^2) * diff(my, psi, 2)
    mzLaplacian = diff(mz, rho, 2) + (1/rho) * diff(mz, rho) + (1/rho^2) * diff(mz, psi, 2)

    # (∇m)^2 term
    quadratic_term = -mx*mxLaplacian - my*myLaplacian - mz*mzLaplacian
    quadratic_term = quadratic_term.simplify_full()

    # (∇^2 m)^2 term
    mLaplacian   = vector([mxLaplacian, myLaplacian, mzLaplacian])
    quartic_term = mLaplacian.dot_product(mLaplacian)
    quartic_term = quartic_term.simplify_full()

    R_values         = np.arange(0.5, 10, 0.1)            # Array of radii used for fitting
    dtheta_dR        = diff(theta, R)                     # dθ/dR
    dtheta_drho      = diff(theta, rho)                   # dθ/dρ
    d2theta_dRdrho   = diff(dtheta_drho, R)               # d^2θ/dρ^2
    dE_2Integrand_dR = 0.5 * diff(rho*quadratic_term, R)  # Derivative of (∇m)^2 term with respect to R
    dE_4Integrand_dR = 0.5 * diff(rho*quartic_term, R)    # Derivative of (∇^2 m)^2 term with respect to R

    dw_width    = 1.4  # Domain wall width
    rho_cutoff  = 50   # Cutoff in integration over radius rho for fitting

    times = np.linspace(0, 1000, 10000)
    np.save(over_directory + '/General/times', times)


    ##################
    # Fit Quantities #
    ##################

    print('Fitting quantities')

    # G12 is the G_{Rη} from the manuscript
    # Gamma11 and Gamma22 are as in the manuscript where 1 is R and 2 is η
    # F_rex is -dU_ex / dR

    # What I call Xi are parameters that show up in the analytical expressions for the helicity-dependent
    # parts of the generalised forces for the skyrmion profile ansatz we use
    
    # Xi1 is the integral of cos(2θ) dθ / dR over ρ
    # Xi2 is the integral of ρ * d^2θ / dRdρ over ρ
    # Xi3 is the integral of cos(θ)sin(θ) over ρ
    # Xi4 is the integral of ρ * dθ / dρ over ρ

    G12 = np.zeros_like(R_values, dtype=float)
    Gamma11 = np.zeros_like(R_values, dtype=float)
    Gamma22 = np.zeros_like(R_values, dtype=float)
    F_Rex = np.zeros_like(R_values, dtype=float)
    Xi1 = np.zeros_like(R_values, dtype=float)
    Xi2 = np.zeros_like(R_values, dtype=float)
    Xi3 = np.zeros_like(R_values, dtype=float)
    Xi4 = np.zeros_like(R_values, dtype=float)

    for i in range(len(R_values)):
        G12[i] = numerical_integral(rho * np.sin(theta(R=R_values[i], w=dw_width, m=1)) * dtheta_dR(R=R_values[i], w=dw_width, m=1), 0, rho_cutoff)[0]
        Gamma11[i] = numerical_integral(rho * dtheta_dR(R=R_values[i], w=dw_width, m=1)^2, 0, rho_cutoff)[0]
        Gamma22[i] = numerical_integral(rho * sin(theta(R=R_values[i], w=dw_width, m=1))^2, 0, rho_cutoff)[0]
        F_Rex[i] = numerical_integral(-dE_4Integrand_dR(R=R_values[i], w=dw_width, m=1) + dE_2Integrand_dR(R=R_values[i], w=dw_width, m=1), 0, rho_cutoff)[0]
        Xi1[i] = numerical_integral(cos(2*theta(R=R_values[i], w=dw_width, m=1)) * dtheta_dR(R=R_values[i], w=dw_width, m=1), 0, rho_cutoff)[0]
        Xi2[i] = numerical_integral(rho * d2theta_dRdrho(R=R_values[i], w=dw_width, m=1), 0, rho_cutoff)[0]
        Xi3[i] = numerical_integral(cos(theta(R=R_values[i], w=dw_width, m=1)) * sin(theta(R=R_values[i], w=dw_width, m=1)), 0, rho_cutoff)[0]
        Xi4[i] = numerical_integral(rho * dtheta_drho(R=R_values[i], w=dw_width, m=1), 0, rho_cutoff)[0]
        

    G12FitParams = curve_fit(linear_fit, R_values, G12, p0=[2., 1.])[0]
    Gamma11FitParams = curve_fit(inverse_linear_fit, R_values, Gamma11, p0=[1.17, 1.5, 0, 0])[0]
    Gamma22FitParams = curve_fit(linear_fit, R_values, Gamma22, p0=[3., -0.5])[0]
    F_RexFitParams = curve_fit(inverse_fourth_order_fit, R_values, F_Rex, p0=[1.17, 1.5, 1, 1, 1, 1])[0]
    Xi1FitParams = curve_fit(inverse_fourth_order_fit, R_values, Xi1, p0=[1.17, 1.5, 1, 1, 1, 1])[0]
    Xi2FitParams = curve_fit(inverse_fourth_order_fit, R_values, Xi2, p0=[1.17, 1.5, 1, 1, 1, 1])[0]
    Xi3FitParams = curve_fit(inverse_fourth_order_fit, R_values, Xi3, p0=[1.17, 1.5, 1, 1, 1, 1])[0]
    Xi4FitParams = curve_fit(linear_fit, R_values, Xi4, p0=[1.17, 1.5])[0]


    ########################
    # Fit Energy Integrals #
    ########################

    # Fit integrals of energy density over space to obtain fitting parameters to use during time integration of collective coordinates

    print('Fitting energies')

    E_exchange = np.zeros_like(R_values, dtype=float)
    E_magnetic_integral = np.zeros_like(R_values, dtype=float)
    E_electric_integral = Xi3 + Xi4

    for i in range(len(R_values)):
        E_exchange[i]          = numerical_integral(-0.5*rho*quadratic_term(R=R_values[i], w=dw_width, m=1) + 0.5*rho*quartic_term(R=R_values[i], w=dw_width, m=1), 0, rho_cutoff)[0]
        E_magnetic_integral[i] = numerical_integral(rho*(mz(R=R_values[i], w=dw_width, m=1)-1), 0, rho_cutoff)[0]

    E_exchangeFitParams          = curve_fit(inverse_quadratic_linear_fit, R_values, E_exchange, p0=[1.17, 1.5, 0, 1])[0]
    E_magnetic_integralFitParams = curve_fit(quadratic_fit, R_values, E_magnetic_integral, p0=[1.17, 1.5, 0])[0]
    E_electric_integralFitParams = curve_fit(linear_fit, R_values, E_electric_integral, p0=[1.17, 1.5])[0]

    # Array of electric field amplitude and angular frequency, each pair of which the collective coordinates are integrated
    E0_values = np.arange(0., 2.05, 0.05)
    omega_values = np.arange(0., 4.02, 0.02)

    Bz = 1.       # Magnetic field
    alpha = 0.01  # Gilbert damping constant

    # Obtain the radius that minimises the energy with zero electric field
    R_values_finer = np.linspace(np.min(R_values), np.max(R_values), 1000)
    E_values_finer = inverse_quadratic_linear_fit(R_values_finer, E_exchangeFitParams[0], E_exchangeFitParams[1], E_exchangeFitParams[2], E_exchangeFitParams[3]) \
        -Bz * quadratic_fit(R_values_finer, E_magnetic_integralFitParams[0], E_magnetic_integralFitParams[1], E_magnetic_integralFitParams[2])
    optimal_skyrmion_radius = R_values_finer[np.argmin(E_values_finer)]

    # Output parameters such as the calculated fitting constants to a file (as they are needed for some of the data analysis scripts)
    with open(over_directory + '/General/CCParams.txt', 'w+') as paramsfile:
        paramsfile.write('Bz=' + str(Bz) + '\n')
        paramsfile.write('dw_width=' + str(dw_width) + '\n')
        paramsfile.write('E_exchangeFitParams=' + ','.join(E_exchangeFitParams.astype(str)) + '\n')
        paramsfile.write('E_magnetic_integralFitParams=' + ','.join(E_magnetic_integralFitParams.astype(str)) + '\n')
        paramsfile.write('E_electric_integralFitParams=' + ','.join(E_electric_integralFitParams.astype(str)) + '\n')
        paramsfile.write('alpha=' + str(alpha) + '\n')
        paramsfile.write('F_RexFitParams=' + ','.join(F_RexFitParams.astype(str)) + '\n')
        paramsfile.write('Gamma11FitParams=' + ','.join(Gamma11FitParams.astype(str)) + '\n')
        paramsfile.write('Gamma22FitParams=' + ','.join(Gamma22FitParams.astype(str)) + '\n')
        paramsfile.write('G12FitParams=' + ','.join(G12FitParams.astype(str)) + '\n')
        paramsfile.write('Xi1FitParams=' + ','.join(Xi1FitParams.astype(str)) + '\n')
        paramsfile.write('Xi2FitParams=' + ','.join(Xi2FitParams.astype(str)) + '\n')
        paramsfile.write('Xi3FitParams=' + ','.join(Xi3FitParams.astype(str)) + '\n')
        paramsfile.write('Xi4FitParams=' + ','.join(Xi4FitParams.astype(str)) + '\n')
        paramsfile.write('OptimalSkyrmionRadius=' + str(optimal_skyrmion_radius) + '\n')


    print('Integrating collective coordinates')

    # Perform the collective coordinate integrations in parallel
    with Pool(multiprocessing.cpu_count(), initialiser, initargs=(over_directory, times, omega_values, Bz, E0_values, alpha, R_values,
        E_exchangeFitParams, E_magnetic_integralFitParams, E_electric_integralFitParams,
        Gamma11FitParams, Gamma22FitParams, G12FitParams, F_RexFitParams,
        Xi1FitParams, Xi2FitParams, Xi3FitParams, Xi4FitParams)) as p:

        p.map(evolve, [Eidx for Eidx in range(len(E0_values))])
