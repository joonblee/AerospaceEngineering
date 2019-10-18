#!/bin/bash
# source respective.sh $1

if [ -z "$1" ]; then
  echo 'If you want to run this code, try below:'
  echo 'source respective.sh $1'
  echo '$1: mass sign of 2nd imperfection'
  exit
else
  touch temp

  MASS=0.1
  PHI=0.5*pi

  echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' >> temp
  echo '%% options                      %%' >> temp
  echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' >> temp
  echo '' >> temp
  if [ -z "$3" ]; then
    echo "m_i     = [$2${MASS}];" >> temp
    echo "phi_i   = [${PHI}];"    >> temp
    echo "theta_i = [$1];"         >> temp

	echo "SIGN_i  = ['$2'];" >> temp
  else
    echo "m_i     = [$2${MASS} $4${MASS}];" >> temp
    echo "phi_i   = [${PHI} ${PHI}];"       >> temp
    echo "theta_i = [$1 $3];"                >> temp

	echo "SIGN_i  = ['$2' '$4'];" >> temp
  fi

  echo '' >> temp

  cat temp main.m >> run.m
  rm temp

  echo "Run code"
  echo "matlab run"

  matlab run > log < /dev/null

  # mv run.m run.m.save
  rm run.m
fi

echo "Finished"

