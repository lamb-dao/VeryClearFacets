#!/bin/bash

#creates sym links in target dir


sou="/home/holdens/holdensQNAP/LIBS/gbs_snp/cassetteGBS/morex/2020/"
ddd="NS.*"
f="variantsMerged.vcf"
tar="$(pwd)/../inputBatch"

mkdir -p "${tar}"

for d in  ${sou}/${ddd}; do

  echo "Enter: ${d}"
  (cd "${d}" || exit

  echo "Build name for: ${f}"
  n=$(basename ${d} | sed s:NS.*Brandon_::g)
  n+=".vcf"
  echo ${n}

  echo "Create renamed links in ${tar}"
  ln -ifs ${d}/${f} ${tar}/${n}
  )
done
