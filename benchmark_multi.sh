#! /bin/bash

nRuns=10
param2Tune="graspChance"
tuningVars="0.1" #0.05 0.025 0.01"
inputDir="instances"
execPath="./main"
solverFixedArgs="-m nn --graspType almostbest --loglvl log --round"
outDir="runs"
outFname="outFullFname"

declare -a subDirs=(
    "0-80"
    "100-200"
    "220-320"
    "400-500"
    "500-800"
    "1000-1440"
    "1570-2400"
    "3000-6000"
    "7000-20000"
)

declare -a tlim=(
    2
    4
    6
    10
    15
    20
    25
    30
    45
)

for i in ${!subDirs[@]}; do
    #echo "${subDirs[$i]} has timeLimit set to ${tlim[$i]}
    mkdir $outDir/${subDirs[$i]} -p

    solverArgs="$solverFixedArgs -t ${tlim[$i]}"
    outTemplateFname="$outDir/${subDirs[$i]}/$outFname${tlim[$i]}sec"
    outCostsFname="$outTemplateFname-costs.csv"
    outRuntimesFname="$outTemplateFname-runtimes.csv"
    outIterCountFname="$outTemplateFname-iterCount.csv"
    fullInputPath="$inputDir/${subDirs[$i]}"
    
    echo benchmark.py --execPath $execPath -n $nRuns --solverExtraArgs "$solverArgs" --inputDir $fullInputPath --param2Tune $param2Tune --tuningVars $tuningVars --saveCosts $outCostsFname --saveRuntimes $outRuntimesFname --saveIterCount $outIterCountFname
    python benchmark.py --execPath $execPath -n $nRuns --solverExtraArgs "$solverArgs" --inputDir $fullInputPath --param2Tune $param2Tune --tuningVars $tuningVars --saveCosts $outCostsFname --saveRuntimes $outRuntimesFname --saveIterCount $outIterCountFname
done
