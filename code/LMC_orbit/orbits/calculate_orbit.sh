#!/bin/sh
# Script to compute the orbit of the LMC
snap_name=MWfLMC3Hdf1
N_i=0
N_f=24
path=../../../data/LMCMW/MWmLMC3/Hernquist/df1/
out_name=orb_Herndf1.txt
N_MW=1000000
N_LMC=500000

python ../../CM/orbit_cm.py $snap_name $N_i $N_f $path $out_name $N_MW $N_LMC
