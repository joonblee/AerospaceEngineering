#!/bin/bash
# source run.sh 1

if [[ ($1 != '1') && ($1 != '2') ]]; then
  echo "[ERROR] \$1 should be '1 or 2'."
  echo "Stop running code."
  # exit 1: Close the terminal
  return 0
fi

DIR=figures_Project$1

if [ -d "$DIR" ]; then # -f: DIRs
  echo "'$DIR' already exists. Create '$DIR.backup'."
  rm "$DIR.backup"
  mv "$DIR" "$DIR.backup"
else
  echo "There is no $DIR directory."
fi

mkdir $DIR
echo "'$DIR' is created."

matlab Project1 > log.Project1.m < /dev/null 

echo "Finished"

