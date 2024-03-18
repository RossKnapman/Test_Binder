"""
Title: Spacetime magnetic hopfions: from internal excitations and braiding of skyrmions

Authors and Affiliations: Ross Knapman[1,2,3], Timon Tausendpfund[1], Sebastián A. Díaz[2,4], Karin Everschor-Sitte[2,3]
    1 Institute of Physics, Johannes Gutenberg University Mainz, 55128 Mainz, Germany
    2 Faculty of Physics, University of Duisburg-Essen, 47057 Duisburg, Germany
    3 Center for Nanointegration Duisburg-Essen (CENIDE), University of Duisburg-Essen, 47057 Duisburg, Germany
    4 Department of Physics, University of Konstanz, 78457 Konstanz, Germany

Purpose: The various fits used to speed up collective coordinate integration and energy calculation.

Written 2023 by Ross Knapman <ross.knapman@uni-due.de>
"""

def linear_fit(R, m, c):
    return m*R + c

def inverse_linear_fit(R, a, b, c, d):
    return a/(R-c) + b*(R-c) + d

def inverse_quadratic_linear_fit(R, a, b, c, d):
    return a/R**2 + b/R + c + d*R

def quadratic_fit(R, a, b, c):
    return a*R**2 + b*R + c

def inverse_fourth_order_fit(R, a, b, c, d, e, f):
    return a/R**4 + b/R**3 + c/R**2 + d/R + e*R + f
