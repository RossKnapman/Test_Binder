"""
Title: Spacetime magnetic hopfions: from internal excitations and braiding of skyrmions

Authors and Affiliations: Ross Knapman[1,2,3], Timon Tausendpfund[1], Sebastián A. Díaz[2,4], Karin Everschor-Sitte[2,3]
    1 Institute of Physics, Johannes Gutenberg University Mainz, 55128 Mainz, Germany
    2 Faculty of Physics, University of Duisburg-Essen, 47057 Duisburg, Germany
    3 Center for Nanointegration Duisburg-Essen (CENIDE), University of Duisburg-Essen, 47057 Duisburg, Germany
    4 Department of Physics, University of Konstanz, 78457 Konstanz, Germany

Purpose: Read specified parameters given a parameters file (such as Micromagnetics/SimulationsAndAnalysis/MicromagneticParams).

Written 2023 by Ross Knapman <ross.knapman@uni-due.de>
"""

def read_quantity(quantity, paramsfile):

    """
    Takes in a quantity key for the parameters and returns the value.

    Args:
        quantity (str): The quantity desired, e.g. I1 (which is first order exchange stiffness)
        paramsfiles (str): Path to the file `MicromagneticParams`, in which these data are stored

    Returns:
        The corresponding quantity read from the file

    """

    with open(paramsfile) as f:
        for line in f.readlines():
            if line.split('=')[0] == quantity:
                return float(line.split('=')[1])
    
    raise ValueError('Invalid quantity key given:', quantity)
