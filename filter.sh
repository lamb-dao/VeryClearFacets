#!/usr/bin/env bash
#strict mode
set -euo pipefail
IFS=$'\n\t'

# snp filter with multiplot support via R and gnuplot
echo "///=======|||============="
echo "//========||=============="
echo "/=========|==============="
echo " BEGIN: VCF(very clear facets)filter V0.2.5"
echo "========---==============/"
echo "=========--=============//"
echo "==========-============///"

# Requires input dir containing one uncompressed merged.vcf file
# as processed by cassetteGBS
# requires checkpoints dir
# Requires an environment with
# up to date R version
# vcftools
# bcftools
# vcflib=1.0.3
# python = 3.11
# tabixpp=1.1.0
# pdflatex

## set flag to skip interactive visualizations
runMode=all
flag="${1:-default}"
if [ "${flag}" == "batch" ]; then
    runMode="${1}"
fi

## Setup locations
home=$(pwd)
checkpoints=${home}/checkpoints
input=${home}/input
data=${home}/data
scripts=${home}/scripts
output=${home}/output
# TODO check for user folder structure is correct

#### functions
function doStep () {
    # doStep stepName, to check if a step is incomplete, returns true if step is required
    ( cd ${checkpoints}
      n=$(ls -1 $1 2>/dev/null | wc -l)
      if [ "$n" -ne 1 ]; then
          #do it
          echo "check: $1 requires completion"
          return 0
      else
          #it is done
          echo "check: $1 is marked complete"
          return 1
    fi
    )
}

function userApprove () {
    #batch mode, no interactive user,  assume approval for all step
    if [ "${runMode}" == "batch" ]; then
        echo "Continuing"
        return 0
    fi
    # user continue or exit
    read -p "To continue enter y: " continue
    if [ "${continue}" != "y" ]; then
        echo "Exiting"
        exit 0
    fi
}

function markStep () {
    # markStep add stepname, to mark complete
    # markStep del stepname, to remove mark
    ( cd ${checkpoints}
      if [ "${1}" == "add" ]; then
          echo "check: $2 will be marked complete"
          touch $2
      fi
      if [ "${1}" == "del" ]; then
          echo "check: $2 will be UN marked"
          rm -f $2
      fi
    )
}

function compareFiles() {
    #compareFiles size file fileAfter
    #compareFiles line file fileAfter
    #compareFiles vars file fileAfter
    type=$1
    bef=$2
    aft=$3
    #
    if [ "${type}" == "size" ]; then
        befSize=$(stat -c%s "$bef")
        aftSize=$(stat -c%s "$aft")
        hSize=$(ls -alh -Bg ${aft} | awk '{print $4}')
        #
        difSize=$(( befSize-aftSize ))
        difPerc=$( echo "scale=2; ${difSize}*100/${befSize}" | bc -l )
        echo "The file size was reduced by ${difPerc} %, ${aft} size is: ${hSize}"
    fi
    # returns an int of file reduction
    if [ "${type}" = "vars" ]; then
        befSize=$(vcfstats ${bef} | grep '^snps:' | awk '{print $2}')
        aftSize=$(vcfstats ${aft} | grep '^snps:' | awk '{print $2}')
        #
        difSize=$(( befSize-aftSize ))
        difPerc=$( echo "scale=9; ${difSize}*100/${befSize}" | bc -l )
        vars="$(echo "scale=0; ${difPerc}*10/10" | bc -l)"
        echo $vars
    fi
    # if arg 1 is "line" do the long lines calculation
    if [ "${type}" == "line" ]; then
        befLine=$(cat $bef | wc -l)
        aftLine=$(cat $aft | wc -l)
        difLine=$(( befLine-aftLine ))
        difPerc=$( echo "scale=2; ${difLine}*100/${befLine}" | bc -l )
        nVar=$( echo "scale=0; ${aftLine}/2" | bc -l )
        echo "Line count reduced by ${difPerc} %"
    fi
}

function varCount () {
    # varCount fileafter
    aft=$1
    aftSize=$(vcfstats $aft | grep '^snps:' | awk '{print $2}')
    echo "${aftSize}"
}


function multiPlot () {
    #batch mode, no interactive user, skip visualization
    if [ "${runMode}" == "batch" ]; then
        echo "Continuing"
        return 0
    fi
    #usage multiPlot elementCode elementName file.vcf
    elementCode=$1
    elementName=$2
    file=$3
    ../scripts/multiPlot.R ${elementCode} ${elementName} ${file}
}

function crossPlot () {
    #batch mode, no interactive user, skip visualization
    if [ "${runMode}" == "batch" ]; then
        echo "Continuing"
        return 0
    fi
    #crossPlot step element in.vcf out.vcf
    step=$1
    elem="${2}"
    bef="$3"
    aft="$4"
    #drop filter from header
    cat "${bef}" | grep -v '##FILTER' > XP_bef.vcf
    cat "${aft}" | grep -v '##FILTER' > XP_aft.vcf
    bef="XP_bef.vcf"
    aft="XP_aft.vcf"
    #extract element from header
    elemSearch="=${elem},"
    headLine=$(head -500 ${bef} \
                   | { grep '#' || test $? = 1; } \
                   | { grep "${elemSearch}" || test $? = 1; })
    #test no match at all
    wordN=$(echo "${headLine}" | wc -w)
    if [ "${wordN}" -eq 0 ]; then
        echo "ERROR: No matching element found in VCF"
        exit 1
    fi
    #test excessive matches
    lineN=$(echo "${headLine}" | wc -l)
    if [ "${lineN}" -gt 1 ]; then
        echo "ERROR: Ambiguous element or VCF header, excess matching elements found in VCF"
        echo "${headLine}"
        exit 1
    fi
    #known 1 line matches
    #test type
    formp=$(echo ${headLine} | { grep 'FORMAT=' || test $? = 1; } | wc -w)
    infop=$(echo ${headLine} | { grep 'INFO=' || test $? = 1; } | wc -w)
    #test if in INFO or FORMAT, build query string
    if [ "${infop}" -eq 0 ]; then
        # must be format
        qStr="[%POS\t%${elem}\n]"
    elif [ "${formp}" -eq 0 ]; then
        # must be info
        qStr="%POS\t%${elem}\n"
    else
        #cant be both, one must be 0
        echo "ERROR: in crossPlot an undefined header has been encountered"
        exit 1
    fi
    # make data files
    bcftools query -f "${qStr}" ${bef} -o XP_befData.txt
    bcftools query -f "${qStr}" ${aft} -o XP_aftData.txt
    # drop lines if value 0
    awk '{if ($2!=0) print $0}' XP_befData.txt > XP_befData_n0.txt
    awk '{if ($2!=0) print $0}' XP_aftData.txt > XP_aftData_n0.txt
    #create script
    gnuPlotScript="XP_plot.gp"
    cat << EOF > $gnuPlotScript
    set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400
    set output "XP_${step}_histo.png"
    set style fill transparent solid 0.7 noborder
    set grid nopolar
    set grid noxtics nomxtics ytics nomytics noztics nomztics nortics nomrtics \
    nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
    set grid layerdefault   lt 0 linecolor 0 linewidth 0.500,  lt 0 linecolor 0 linewidth 0.500
    set style data lines
    #set xtics border in scale 1,0.5 nomirror norotate autojustify
    #set xtics  norangelimit 1
    set title "${step}"
    set xlabel "${elem}"
    set ylabel "Frequency"
    set xrange [ * : * ] noreverse writeback
    set x2range [ * : * ] noreverse writeback
    set yrange [ 0.00000 : * ] noreverse nowriteback
    set y2range [ * : * ] noreverse writeback
    set zrange [ * : * ] noreverse writeback
    set cbrange [ * : * ] noreverse writeback
    set rrange [ * : * ] noreverse writeback
    set jitter overlap 1  spread 0.5  wrap 0
    set colorbox vertical origin screen 0.9, 0.2 size screen 0.05, 0.6 front  noinvert bdefault
    NO_ANIMATION = 1
    plot "XP_befData_n0.txt" using 2 bins=100 with boxes lc rgb "dark-violet" title 'Before (100 bins)', \
         '' using 2:(1) smooth kdensity bandwidth .5 lw 2 lc rgb "dark-violet" title 'Before (KDE bandwidth:0.5)', \
         "XP_aftData_n0.txt" using 2 bins=100 with boxes lc rgb "forest-green" title 'After (100 bins)', \
         '' using 2:(1) smooth kdensity bandwidth .5 lw 2 lc rgb "forest-green" title 'After (KDE bandwidth:0.5)'
EOF
    #run script
    gnuplot $gnuPlotScript
}

#===================
# clean data and output for this run, leave input and checkpoints
step="sweep"
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    rm -r ${data}
    mkdir -p ${data}
    touch ${data}/.gitkeep
    markStep add $step
fi
#===================

#===================
#Load a single vcf from input to data
step="copyIP"
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    ( cd ${input}
      n=$(ls -1 *.vcf 2>/dev/null | wc -l)
      if [ "$n" -ne 1 ]; then
          echo "ERROR: there can only be one .vcf file in input"
          exit 1
      else
          vcf="$(pwd)/$(ls -1 *.vcf)"
          # full copy
          # rsync -vah --progress $vcf ${data}/in.vcf
          # link copy
          ln -s $vcf ${data}/in.vcf
          markStep add $step
      fi
    )
fi
#===================

#===================
## subset file for handling
step="subsetFile"
in="in.vcf"
out="${step}.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    #optional only if interactive mode is on, in batch mode use full file
    if [ "${runMode}" == "all" ]; then
        read -p "To filter the ENTIRE variant file enter y, otherwise create a subset: " continue
        if [ "${continue}" = "y" ]; then
            echo "Continue with full variant file"
        else
            samplingRate="0.05"
            echo "Subsampling the variant file for filter testing"
            echo "samplingRate=${samplingRate}"
            ( cd ${data}
              check=$(ls -1 ${in})
              vcfrandomsample -p 777 -r ${samplingRate} ${in} > ${out}
            )
        fi
    fi
    markStep add $step
fi
#===================

#===================
# user can skip subset step, this will handle the cases
step="handleSubset"
inFull="in.vcf"
inSub="subsetFile.vcf"
out="filterStart.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    ( cd ${data}
      if [ -e "${inSub}" ]; then
          echo "User created a subsetfile, this is for quickly testing filter choices"
          ln -s "${inSub}" "${out}"
      elif [ -e "${inFull}" ]; then
          echo "User did not create a subsetfile. Continue the pipeline to filter all variants, this will take time but give the best results."
          ln -s "${inFull}" "${out}"
      else
          echo "ERROR: either ${inFull} or ${inSub} must exist"
          exit 1
      fi
      #materialize the symlink
      rsync -Lvh --progress ${out} ${out}.full
      rm ${out}
      mv ${out}.full ${out}
    )
    markStep add $step
fi
#===================

#===================
# manipulates dataField names which differ due to caller
step="renameFields"
in="filterStart.vcf"
#markStep del $step #force do
#markStep add $step #force skip, it is better to try to get bcf filters to handle high precision filters where names are inconsistent with vcftools filter assumptions. This may not be possible
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    # rename NR to DP
    ( cd ${data}
      #FORMAT=<ID=NR
      sed -i "s/FORMAT=<ID=NR/FORMAT=<ID=DP/g" ${in}
      #:NR
      sed -i "s/:NR/:DP/g" ${in}
      sed -i "s/:NR:/:DP:/g" ${in}
    )
    markStep add $step
fi
#===================

#===================
# retains bi-allelic snps only
step="removeIndels"
in="filterStart.vcf"
out="${step}.recode.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    #
    ( cd ${data}
      vcftools \
          --vcf ${in} \
          --remove-indels \
          --min-alleles 2 \
          --max-alleles 2 \
          --recode \
          --recode-INFO-all \
          --out ${step}
      echo -e "SUMMARY\tSTEP\tPERCENT_CUT\tVARS_KEPT"
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================

#===================
# remove sites below MQx phred scale
step="minQuality"
in="${out}" #grab last steps output
out="${step}.recode.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    #
    ( cd ${data}
      vcftools \
          --vcf ${in} \
          --minQ 30 \
          --recode \
          --recode-INFO-all \
          --out ${step}
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================

#===================
# set to missing any genotypes below GQx (phred score)
step="minQualityGT"
in="${out}" #grab last steps output
out="${step}.recode.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    #
    ( cd ${data}
      vcftools \
          --vcf ${in} \
          --minGQ 30 \
          --recode \
          --recode-INFO-all \
          --out ${step}
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================

#===================
# exclude sites with less than x(integer count) total observations, used to remove noise and erronious sequencing
step="minorAlleleCount"
in="${out}" #grab last steps output
out="${step}.recode.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    #
    ( cd ${data}
      vcftools \
          --vcf ${in} \
          --mac 2 \
          --recode \
          --recode-INFO-all \
          --out ${step}
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================

#===================
# set to missing any genotypes with less than x(integer count) reads
step="minDepth"
in="${out}" #grab last steps output
out="${step}.recode.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    #
    ( cd ${data}
      vcftools \
          --vcf ${in} \
          --minDP 5 \
          --recode \
          --recode-INFO-all \
          --out ${step}
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================

#===================
# set to missing any genotypes with more than x(integer count) reads
step="maxDepth"
in="${out}" #grab last steps output
out="${step}.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    #
    ( cd ${data}
      bcftools filter \
               --threads 24 \
               --e 'DP > 400' \
               --set-GTs . \
               -O v \
               -o "${step}.vcf" \
               ${in}
      rm -f XP_*
      crossPlot ${step} DP ${in} ${out}
      multiPlot DP Depth ${out}
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================

#===================
step="minQualDepth"
in="${out}" #grab last steps output
out="${step}.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    #
    ( cd ${data}
      bcftools filter \
               --threads 24 \
               --e 'QD < 3' \
               -O v \
               -o "${step}.vcf" \
               ${in}
      rm -f XP_*
      crossPlot ${step} QD ${in} ${out}
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================


#===================
# exclude individuals sequenced at less than x(ratio) sites
step="maxMissingBySample"
in="${out}" #grab last steps output
out="${step}.recode.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    ( cd ${data}
      #report on all samples
      vcftools \
          --vcf "${in}" \
          --missing-indv
      #find samples above cutoff
      cutoff=0.88
      awk -v co="${cutoff}" '{if ($5<co) print $1}' out.imiss > out.imiss.keep
      awk -v co="${cutoff}" '{if ($5>=co) print $1}' out.imiss > out.imiss.drop
      #remove samples above cutoff
      vcftools \
                  --vcf "${in}" \
                  --keep out.imiss.keep \
                  --recode \
                  --recode-INFO-all \
                  --out "${step}"
      #report
      echo "DROPPING"
      awk  '{print $1}' out.imiss.drop
      echo -e "LINE\tMISSING"
      sort -t$'\t' -k 5 -n -r out.imiss > out.imiss.sort
      awk  '{print $1"\t"$5}' out.imiss.sort
      echo -e "SUMMARY\t${step}\t$(compareFiles line out.imiss out.imiss.keep)\tkeep $(cat out.imiss.keep | wc -l) individuals"
    )
    markStep add $step
fi
#===================

#===================
# exclude sites sequenced in less than x(ratio) of the individuals
step="maxMissingBySite"
in="${out}" #grab last steps output
out="${step}.recode.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    #
    ( cd ${data}
      vcftools \
              --vcf ${in} \
              --max-missing 0.90 \
              --recode \
              --recode-INFO-all \
              --out ${step}
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================

#===================
step="minorAlleleFreq"
in="${out}" #grab last steps output
out="${step}.recode.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    ( cd ${data}
      vcftools \
          --vcf ${in} \
          --maf 0.05 \
          --recode \
          --recode-INFO-all \
        --out ${step}
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================

#===================
step="allelicBalance"
in="${out}" #grab last steps output
out="${step}.recode.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    ( cd ${data}
      vcftools \
          --vcf ${in} \
          --freq2
      awk '{if ($5>0.25 && $5<0.75) print $1"\t"$2}' out.frq > out.frq.keep
      vcftools \
          --vcf "${in}" \
          --positions out.frq.keep \
          --recode \
          --recode-INFO-all \
          --out "${step}"
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================

#===================
# removes sites that exhibit strandbias
# fishers strand bias test gives probability of strand bias
# converted from p value to phred score, more phred more chance of bias
# pvalue phred
# 1 0    ie.noChanceOfBias
# 0.1 10  10% chance that there is no statistically convincing evidence of bias
step="strandbias"
in="${out}" #grab last steps output
out="${step}.recode.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    #
    ( cd ${data}
      # bcftools filter \
      #          --threads 24 \
      #          --e 'FS > 60' \
      #          -O v \
      #          -o "${out}" \
      #          ${in}
      # in our case we filter by flag as FS field not present
      vcftools \
          --vcf "${in}" \
          --remove-filtered "strandBias" \
          --recode \
          --recode-INFO-all \
          --out "${step}"
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================

#===================
# final filter thin file to target snp count
step="thin"
in="${out}"
out="${step}.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    #samplingRates
    x1="0.01"
    x2="0.1"
    y3=20000 #the target variable count in final file
    ( cd ${data}
      check=$(ls -1 ${in})
      vcfrandomsample -p 777 -r "${x1}" ${in} > ${out}
      y1=$(varCount ${out})
      vcfrandomsample -p 777 -r "${x2}" ${in} > ${out}
      y2=$(varCount ${out})
      #calculate slope, m=y2-y1/x2-x1
      m=$( echo "scale=9; ((${y2})-(${y1}))/((${x2})-(${x1}))" | bc -l )
      #calculate intercept, b=y2-m*x2
      b=$( echo "scale=9; (${y2})-((${m})*(${x2}))" | bc -l )
      #solve for sampling rate, x3=y3-b/m
      x3=$( echo "scale=5; ((${y3})-(${b}))/(${m})" | bc -l )
      #calculate final thinned file
      vcfrandomsample -p 777 -r "${x3}" ${in} > ${out}
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================

#===================
# extract summary lines for inter-run comparison
step="summary"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    cat log*.txt | grep SUMMARY > ${data}/summary.raw.txt
    ( cd ${data}
      tac summary.raw.txt \
          | awk '!seen[$2]++' \
          | tac - \
          > summary.txt
    )
    markStep add $step
fi
#===================

#===================
step="finalize"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    cp log*.txt ${data}
    ( cd ${data}
      #include real file name for tracking
      runName=$(cd ${input}; n="*.vcf"; echo ${n} | sed s/.vcf//)
      op="filterOP_${runName}"
      #make unique dir
      rm -r ${op} || echo ""
      mkdir -p ${op}/logs
      #package outputs
      cp summary.txt ${op}
      cp thin.vcf ${op}
      cp strandbias.recode.vcf ${op}/full.vcf
      cp log*.txt ${op}/logs
      # cp *.log ${op}/logs #error here
      # cp *.pdf ${op}
      # cp *.png ${op}
      #push outputs
      mv -f ${op} ${output}
    )
    markStep add $step
fi
#===================

#===================
step="sweepClose"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    ( cd ${input}
      rm *.vcf || echo ""
      cd ${checkpoints}
      rm * || echo ""
      cd ${data}
      rm -r * || echo ""
      cd ..
      rm log*.txt
    )
fi
#===================

exit 0
###################################
#TODO
finetune each step, use crossplot
note site renaming is further along in cassetteGBS

#===================
step="templateStep"
in="${out}" #grab last steps output
out="${step}.vcf"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    ( cd ${data}
      echo "...step action"
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
      rm -f XP_*
      crossPlot ${step} ELEM... ${in} ${out}
      rm -f MP_*
      multiPlot ELEM... ${step} ${out}
    )
    markStep add $step
fi
#===================

#======================================
#filter cycle...
make fully generic
filterStep="..."
in="filterStart.vcf"
out="${filterStep}.recode.vcf"
#===================
step="${filterStep}VisPre"
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    ( cd ${data}
      rm -f MP_*
      multiPlot DP ${step} ${in}
    )
    markStep add $step
fi
#===================
step="${filterStep}"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    ( cd ${data}
      vcftools \
          --vcf ${in} \
          ...someStep...
          --...?remove-filtered-all \ #want just this step
          --recode \
          --recode-INFO-all \
          --out ${step}
      compareFiles size ${in} ${out}
      compareFiles line ${in} ${out}
      echo -e "SUMMARY\t${step}\t$(compareFiles vars ${in} ${out})\t$(varCount ${out})"
    )
    markStep add $step
fi
#===================
step="${filterStep}VisPost"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "BEGIN: ${step}"
    userApprove
    ( cd ${data}
      multiPlot DP ${step} ${out}
    )
    markStep add $step
fi
#===================
step="${filterStep}Finalize"
#markStep del $step #force do
#markStep add $step #force skip
if( doStep $step ); then
    echo "Examine post filter visualization before continuing."
    userApprove
    echo "BEGIN: ${step}"
    echo -e "If the results of the filter are acceptable then continue to the next step,\n\totherwise exit now. This filter will run again on the next call.\n\tThe filter settings you may want to modify are before line number: ${LINENO}\n\tOnly enter y to finalize this filter"
    echo "To continue enter y: "
    read continue
    if [ "${continue}" != "y" ]; then
        ( cd ${checkpoints}
          # delete checkpoints after previs
          rm ${filterStep} ${filterStep}VisPost
          # remove all old MP*VisPost* files
          cd ${data}
          rm MP*VisPost*
        )
        echo "Exiting, to repeat run again"
        exit 0
    fi
    markStep add $step
fi
#===================
#======================================


## these fields are available in platypus variant files
      nodata=(FS START ReadPosRankSum SC Size)
INFO=<ID=FR,Number=.,Type=Float,Description="Estimated population frequency of variant">
INFO=<ID=MMLQ,Number=1,Type=Float,Description="Median minimum base quality for bases around variant">
INFO=<ID=TCR,Number=1,Type=Integer,Description="Total reverse strand coverage at this locus">
INFO=<ID=HP,Number=1,Type=Integer,Description="Homopolymer run length around variant locus">
INFO=<ID=WE,Number=1,Type=Integer,Description="End position of calling window">
INFO=<ID=Source,Number=.,Type=String,Description="Was this variant suggested by Playtypus, Assembler, or from a VCF?">
#INFO=<ID=FS,Number=.,Type=Float,Description="Fisher's exact test for strand bias (Phred scale)">
INFO=<ID=WS,Number=1,Type=Integer,Description="Starting position of calling window">
INFO=<ID=PP,Number=.,Type=Float,Description="Posterior probability (phred scaled) that this variant segregates">
INFO=<ID=TR,Number=.,Type=Integer,Description="Total number of reads containing this variant">
INFO=<ID=NF,Number=.,Type=Integer,Description="Total number of forward reads containing this variant">
INFO=<ID=TCF,Number=1,Type=Integer,Description="Total forward strand coverage at this locus">
INFO=<ID=NR,Number=.,Type=Integer,Description="Total number of reverse reads containing this variant">
INFO=<ID=TC,Number=1,Type=Integer,Description="Total coverage at this locus">
INFO=<ID=END,Number=.,Type=Integer,Description="End position of reference call block">
INFO=<ID=MGOF,Number=.,Type=Integer,Description="Worst goodness-of-fit value reported across all samples">
INFO=<ID=SbPval,Number=.,Type=Float,Description="Binomial P-value for strand bias test">
#INFO=<ID=START,Number=.,Type=Integer,Description="Start position of reference call block">
#INFO=<ID=ReadPosRankSum,Number=.,Type=Float,Description="Mann-Whitney Rank sum test for difference between in positions of variants in reads from ref and alt">
INFO=<ID=MQ,Number=.,Type=Float,Description="Root mean square of mapping qualities of reads at the variant position">
INFO=<ID=QD,Number=1,Type=Float,Description="Variant-quality/read-depth for this variant">
#INFO=<ID=SC,Number=1,Type=String,Description="Genomic sequence 10 bases either side of variant position">
INFO=<ID=BRF,Number=1,Type=Float,Description="Fraction of reads around this variant that failed filters">
INFO=<ID=HapScore,Number=.,Type=Integer,Description="Haplotype score measuring the number of haplotypes the variant is segregating into in a window">
#INFO=<ID=Size,Number=.,Type=Integer,Description="Size of reference call block">

FILTER=<ID=GOF,Description="Variant fails goodness-of-fit test.">
FILTER=<ID=hp10,Description="Flanking sequence contains homopolymer of length 10 or greater">
FILTER=<ID=REFCALL,Description="This line represents a homozygous reference call">
FILTER=<ID=badReads,Description="Variant supported only by reads with low quality bases close to variant position, and not present on both strands.">
FILTER=<ID=alleleBias,Description="Variant frequency is lower than expected for het">
FILTER=<ID=Q20,Description="Variant quality is below 20.">
FILTER=<ID=HapScore,Description="Too many haplotypes are supported by the data in this region.">
FILTER=<ID=MQ,Description="Root-mean-square mapping quality across calling region is low.">
FILTER=<ID=QD,Description="Variants fail quality/depth filter.">
FILTER=<ID=SC,Description="Variants fail sequence-context filter. Surrounding sequence is low-complexity">
FILTER=<ID=QualDepth,Description="Variant quality/Read depth ratio is low.">
FILTER=<ID=strandBias,Description="Variant fails strand-bias filter">

FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">
FORMAT=<ID=GQ,Number=.,Type=Integer,Description="Genotype quality as phred score">
FORMAT=<ID=GOF,Number=.,Type=Float,Description="Goodness of fit value">
>>>> NR > DP
FORMAT=<ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">
FORMAT=<ID=GL,Number=.,Type=Float,Description="Genotype log10-likelihoods for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites">
FORMAT=<ID=NV,Number=.,Type=Integer,Description="Number of reads containing variant in this sample">

#site line
R890096.1	162802	.	A	G	2965	PASS
BRF=0.29;FR=0.3122;HP=1;HapScore=1;MGOF=121;MMLQ=25;MQ=60;NF=871;NR=2;PP=2965;QD=20;SC=AGGGAACAACAACTCATGAGG;SbPval=0.55;Source=Platypus;TC=2550;TCF=2547;TCR=3;TR=873;WE=162810;WS=162792
GT:GL:GOF:GQ:NR:NV
0/0:0,-6.9,-78.19:1:69:22:0
0/0:0,-4.52,-52.39:36:45:45:0
0/0:0,-9.64,-115.29:2:96:32:0
1/1:-33.71,-2.33,0:27:23:15:15
0/0:0,-0.34,-4.39:17:5:6:0 ...


# header record of settings
platypusOptions={'assemblyRegionSize': 1500, 'trimReadFlank': 0, 'assembleBadReads': 1, 'bamFiles': ['data/NS.2038.003.B722---B504.GBS_AAFC-Brandon_2022_P16/chr1H/listBam.txt'], 'minVarDist': 9, 'trimSoftClipped': 1, 'minReads': 2, 'qualBinSize': 1, 'refFile': 'cassettes/cassette_NS.2038.003.B722---B504.GBS_AAFC-Brandon_2022_P16/refgenome/chr1H.fasta', 'maxHaplotypes': 50, 'filterVarsByCoverage': 0, 'maxSize': 1500, 'originalMaxHaplotypes': 50, 'skipDifficultWindows': 0, 'parseNCBI': 0, 'skipRegionsFile': None, 'noCycles': 0, 'trimAdapter': 0, 'minPosterior': 5, 'assembleAll': 1, 'trimOverlapping': 1, 'filterDuplicates': 0, 'abThreshold': 0.01, 'minFlank': 5, 'bufferSize': 100000, 'fileCaching': 0, 'useEMLikelihoods': 0, 'coverageSamplingLevel': 30, 'calculateFlankScore': 0, 'logFileName': 'logs/NS.2038.003.B722---B504.GBS_AAFC-Brandon_2022_P16/platypus/chr1H.log', 'nCPU': 8, 'filterReadsWithUnmappedMates': 0, 'qdThreshold': 10, 'maxVariants': 8, 'scThreshold': 0.95, 'filterReadsWithDistantMates': 1, 'maxReads': 500000000.0, 'badReadsWindow': 11, 'genIndels': 1, 'largeWindows': 0, 'minMapQual': 20, 'maxVarDist': 15, 'maxGOF': 20, 'rlen': 250, 'minGoodQualBases': 5, 'refCallBlockSize': 1000, 'countOnlyExactIndelMatches': 0, 'longHaps': 0, 'HLATyping': 0, 'filterReadPairsWithSmallInserts': 1, 'minBaseQual': 20, 'getVariantsFromBAMs': 1, 'genSNPs': 1, 'assemble': 0, 'assemblerKmerSize': 15, 'minVarFreq': 0.002, 'alignScoreFile': '', 'verbosity': 2, 'sourceFile': None, 'compressReads': 0, 'rmsmqThreshold': 20, 'filteredReadsFrac': 0.7, 'outputRefCalls': 0, 'badReadsThreshold': 10, 'hapScoreThreshold': 15, 'regions': None, 'sbThreshold': 0.01, 'output': 'data/NS.2038.003.B722---B504.GBS_AAFC-Brandon_2022_P16/chr1H/variants.vcf', 'assembleBrokenPairs': 0, 'mergeClusteredVariants': 0, 'maxGenotypes': 1275, 'nInd': 94}
