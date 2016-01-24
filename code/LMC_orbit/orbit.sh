#!/bin/bash
path=/rsgrps/gbeslastudents/nicolas/LMCMW/simulations/LMCMW/LR/MW1LMC4/a1/
file=MW1LMC4a1H3*
mvto=../../data/LMCMW/MW1LMC4/a1
path2snap=../../data/LMCMW/MW1LMC4/a1/
snaps=MW1LMC4a1H3
orbit=LMCMW-H3.txt

#scp jngaravitoc@elgato-login.hpc.arizona.edu:$path$file $mvto
python ../CM/orbit_cm.py  $snaps 0 100 $path2snap $orbit  30000 30000
