#!/bin/sh

set -x -e

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

inputfile=$1
indexname=`basename ${inputfile}`
for lcptab_alg in kasai13n kasai9n plcp5n
do
  ./sa_induced.x --verbose --absolute_suftab --indexname ${indexname}.${lcptab_alg} -l ${lcptab_alg} ${inputfile}
done
for suffix in lcp ll2 ll4 prj isa suf
do
  if test -f ${indexname}.kasai13n.${suffix}
  then
    if test -f ${indexname}.kasai9n.${suffix}
    then
      cmp -s ${indexname}.kasai13n.${suffix} ${indexname}.kasai9n.${suffix}
    fi
    if test -f ${indexname}.plcp5n.${suffix}
    then
      cmp -s ${indexname}.kasai13n.${suffix} ${indexname}.plcp5n.${suffix}
    fi
  else
    if test -f ${indexname}.kasai9n.${suffix}
    then
      if test -f ${indexname}.plcp5.${suffix}
      then
        cmp -s ${indexname}.kasai9n.${suffix} ${indexname}.plcp5n.${suffix}
      fi
    fi
  fi
done
for suffix in lcp ll2 ll4 prj isa suf
do
  for lcptab_alg in kasai13n kasai9n plcp5n
  do
    rm -f ${indexname}.${lcptab_alg}.${suffix}
  done
done
