#!/bin/sh

set -x -e

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

inputfile=$1
# The default protein alphabet for the genometools is
# alphabet::amino_acids_gt_order
# while sk-sain uses
# alphabet::amino_acids. So we provide the smap file
# TransProtAlphabetOrder.txt which is compatible with
# alphabet::amino_acids.

if test "${GTDIR}" = ""
then
  echo "$0: environment variable GTDIR not defined => skip test"
  exit 0
fi
if test -d "${GTDIR}" -ne 0
then
  echo "$0: directory ${GTDIR} does not exist => skip test"
  exit 0
fi

is_protein=$("${GTTL}"/scripts/guess_if_protein_seq.py --print_result "${inputfile}")
if test "${is_protein}" = "1"
then
  file_type_option="-smap TransProtAlphabetOrder.txt"
else
  file_type_option=-dna
fi
sfx=$(mktemp --tmpdir=. TMP.XXXXXX) || exit 1
cmd="${GTDIR}/bin/gt encseq encode ${file_type_option} -indexname ${sfx} ${inputfile}"
${cmd}
if test $? -ne 0
then
  echo "$0: FAILURE: ${cmd}"
  exit 1
fi

sain=$(mktemp --tmpdir=. TMP.XXXXXX) || exit 1
cmd="./sa_induced.x --indexname ${sain} --lcptab plcp5n -a ${inputfile}"
${cmd}
if test $? -ne 0
then
  echo "$0: FAILURE: ${cmd}"
  exit 1
fi

if test -f "${sain}".prj
then
  sizeof_suftab_entry=$(grep '^sizeof_suftab_entry' "${sain}".prj | cut -f 2)
else
  echo "$0: FAILURE: ${sain}.prj does not exist"
  exit 1
fi

if test "${sizeof_suftab_entry}" -eq 32
then
  suffixerator_size_option=-suftabuint
fi

cmd="${GTDIR}/bin/gt suffixerator -ii ${sfx} -lcp -suf ${suffixerator_size_option}"
${cmd}
if test $? -ne 0
then
  echo "$0: FAILURE: ${cmd}"
  exit 1
fi

for suffix in suf lcp
do
  cmd="cmp -s ${sfx}.${suffix} ${sain}.${suffix}"
  ${cmd}
  if test $? -ne 0
  then
    echo "$0: FAILURE: ${cmd}"
    exit 1
  fi
  rm -f "${sfx}".${suffix} "${sain}".${suffix}
done
cmd="intread.x 8 ${sfx}.llv"
${cmd} > "${sfx}".llv.txt
if test $? -ne 0
then
  echo "$0: FAILURE: ${cmd}"
  exit 1
fi
if test -f "${sain}".llv
then
  cmd="intread.x 4 ${sain}.llv"
  ${cmd} > "${sain}".llv.txt
  if test $? -ne 0
  then
    echo "$0: FAILURE: ${cmd}"
    exit 1
  fi
  cmd="diff --strip-trailing-cr ${sfx}.llv.txt ${sain}.llv.txt"
  ${cmd}
  if test $? -ne 0
  then
    echo "$0: FAILURE: ${cmd}"
    exit 1
  fi
fi
for suffix in ll2 ll4 lcp prj llv.txt bsf
do
  rm -f "${sain}".${suffix}
done
rm -f "${sain}"
for suffix in llv des sds md5 esq ssp prj llv.txt
do
  rm -f "${sfx}".${suffix}
done
rm -f "${sfx}"
