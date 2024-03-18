#!/bin/bash

# Title: Spacetime magnetic hopfions: from internal excitations and braiding of skyrmions
#
# Authors and Affiliations: Ross Knapman[1,2,3], Timon Tausendpfund[1], Sebastián A. Díaz[2,4], Karin Everschor-Sitte[2,3]
#    1 Institute of Physics, Johannes Gutenberg University Mainz, 55128 Mainz, Germany
#    2 Faculty of Physics, University of Duisburg-Essen, 47057 Duisburg, Germany
#    3 Center for Nanointegration Duisburg-Essen (CENIDE), University of Duisburg-Essen, 47057 Duisburg, Germany
#    4 Department of Physics, University of Konstanz, 78457 Konstanz, Germany
#
# Purpose: Read in the values of E0 and omega (eletric field amplitude and angular frequency) from SimulationE0Omegas.txt
# and, for each, run a simulation using parameters defined in MicromagneticParams. Once the simulation has been run, analyse the results
# using the Python script AnalyseSimulations.py
#
# Written 2023 by Ross Knapman <ross.knapman@uni-due.de>

output_over_dir='MicromagneticsData'
original_dir=$(pwd -P)

mkdir -p $output_over_dir

# Load config parameters
cp MicromagneticParams $output_over_dir
source MicromagneticParams

while IFS="," read -r E0 omega
do

    output_dir=$output_over_dir/E0${E0}omega$omega
    mkdir -p $output_dir

    # Replace placeholders in Template.mx3 with values from MicromagneticParams
    sed "s/{{ N }}/$N/" Template.mx3 |\
    sed "s/{{ DiscretisationSize }}/$discretisationsize/" |\
    sed "s/{{ TimeDiscretisation }}/$timediscretisation/" |\
    sed "s/{{ Ms }}/$Ms/" |\
    sed "s/{{ I1 }}/$I1/" |\
    sed "s/{{ I2 }}/$I2/" |\
    sed "s/{{ Lbda }}/$lbda/" |\
    sed "s/{{ DimensionlessField }}/$dimensionlessfield/" |\
    sed "s/{{ RelaxTime }}/$relaxtime/" |\
    sed "s/{{ SimulationTime }}/$simulationtime/" |\
    sed "s/{{ Alpha }}/$alpha/" |\
    sed "s/{{ AbsorbingSize }}/$absorbingsize/" |\
    sed "s/{{ EFieldOmega }}/$omega/" |\
    sed "s/{{ DimensionlessEAmplitude }}/$E0/" > $output_dir/Sim.mx3

    cd $output_dir/

    ../../../ModifiedMuMax3/modified_mumax3_binary -f -o Data Sim.mx3
    python ../../AnalyseSimulations.py .
    cd $original_dir

done < SimulationE0Omegas.txt
