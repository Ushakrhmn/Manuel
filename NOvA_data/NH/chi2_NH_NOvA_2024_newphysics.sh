#!/bin/sh
export LD_LIBRARY_PATH=/eos/user/u/urahaman/SWAN_projects/long-baseline/globes/globes-3.2.18/lib:\$LD_LIBRARY_PATH
cd /afs/cern.ch/user/u/urahaman/Manuel/NOvA_data/NH/

./chi2_NH_NOvA_2024_newphysics $1
