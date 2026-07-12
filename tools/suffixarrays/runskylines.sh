#!/bin/sh

if test ! -d ${PACKAGES}/libdivsufsort/lib
then
  exit 0
else
  export DYLD_LIBRARY_PATH=${PACKAGES}/libdivsufsort/lib
fi

# compare sais lite and libdivsufsort for skyline strings.

param=3
while test $param -ne 25
do
  for direction in "" "--reverse"
  do
    echo "skyline ${direction} ${param}"
    TMPFILE=`mktemp --tmpdir=. TMP.XXXXXX` || exit 1
    ./skyline.py ${direction} ${param} > ${TMPFILE}
    if test $? -ne 0
    then
      echo "FAILURE: ./skyline.py ${direction} ${param}"
      exit 1
    fi
    for lcptab in "" "--lcptab kasai13n --check_suftab" "--absolute_suftab --lcptab kasai9n" "--absolute_suftab --lcptab plcp5n"
    do
      ./sa_induced.x -v --plain_input_format ${lcptab} ${suftab} ${TMPFILE}
      if test $? -ne 0
      then
        echo "FAILURE: ./sa_induced.x -v --plain_input_format ${lcptab} ${suftab} ${TMPFILE}"
        exit 1
      fi
    done
    rm -f ${TMPFILE}
    for suffix in suf ll2 ll4 lcp prj
    do
      rm -f ${TMPFILE}.${suffix}
    done
  done
  param=`expr $param + 1`
done
