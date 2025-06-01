#!/bin/sh

set -e -x

# compare sais lite and libdivsufsort for skyline strings.

param=3
export DYLD_LIBRARY_PATH=${PACKAGES}/libdivsufsort/lib
while test $param -ne 25
do
  for direction in "" "--reverse"
  do
    echo "skyline ${direction} ${param}"
    outputfile=skyline${direction}-${param}.txt
    ./skyline.py ${direction} ${param} > ${outputfile}
    for lcptab in "" "--lcptab kasai13n --check_suftab" "--absolute_suftab --lcptab kasai9n" "--absolute_suftab --lcptab plcp5n"
    do
      ./sa_induced.x -v --plain_input_format ${lcptab} ${suftab} ${outputfile}
    done
    rm -f ${outputfile}
    for suffix in suf ll2 ll4 lcp prj
    do
      rm -f ${outputfile}.${suffix}
    done
  done
  param=`expr $param + 1`
done
