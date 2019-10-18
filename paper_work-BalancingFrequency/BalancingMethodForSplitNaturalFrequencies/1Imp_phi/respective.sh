#!/bin/bash
# source respective.sh ${mode_number} ${theta1} ${mass_sign1} ${theta2} ${mass_sign2}

DIR=figures

if [ -d "$DIR" ]; then # -f: DIRs
  echo "Directory '$DIR' already exists."
else
  echo "There is no '$DIR' directory."
  mkdir $DIR
  echo "'$DIR' is created."
fi

if [ -z "$1" ]; then
  echo 'If you want to run this code, try below:'
  echo 'source respective.sh $1'
  echo '$1: mode number'
else
  rm run.m

  touch temp

  echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' >> temp
  echo '%% options                      %%' >> temp
  echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' >> temp
  echo '' >> temp
  echo "n=$1; %% mode number" >> temp
  echo '' >> temp

  # imperfection #
  MASS=0.1

  # ------------ #
  echo "phi0 = 2*pi/3;" >> temp
  #echo "phi0 = pi/2;" >> temp
  echo "THETA = 0;" >> temp
  echo "PHI  = -phi0 + 2*[0:1800]/1800*phi0;"        >> temp
  echo "" >> temp

  echo "phi_i   = [PHI]" >> temp
  echo "theta_i = [THETA]" >> temp
  echo "m_i     = [${MASS}]" >> temp
  echo "SIGN_i  = ['$2'];" >> temp

  echo '' >> temp
  #echo 'phi   = pi/2;' >> temp
  echo '' >> temp


  cat temp main.m >> run.m
  rm temp

  echo "Run code"
  echo "matlab run ($1)"

  matlab run > log < /dev/null

  rm run.m
fi

echo "Finished"

