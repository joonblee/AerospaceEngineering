#!/bin/bash
# source run.sh 

DIR=figures

if [ -d "$DIR" ]; then # -f: DIRs
  echo "'$DIR' already exists."
else
  echo "There is no $DIR directory."
  mkdir $DIR
  echo "'$DIR' is created."
fi

matlab fig3 > log < /dev/null

#if [[ ($1 == '2') && ($1 == '4') ]]; then
#  matlab fig$1 > log.fig$1 < /dev/null
#elif [[ ($1 == '3') ]]; then
#  cat fig3_option.m fig3_skel.m > fig3.m
#  matlab fig3 > log.fig3 < /dev/null
#  rm fig3.m
#else
#  echo "[ERROR] \$1 should be 2, 3, or 4."
#  echo "Try 'source run.sh \$1'."
#  echo "Stop running code."
#  # exit 1: Close the terminal
#  return 0
#fi

echo "Finished"

