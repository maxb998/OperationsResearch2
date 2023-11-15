#! /bin/bash

BenchmarkPy="pyScripts/Benchmark.py"

nRuns=3
param2Tune="graspChance"
tuningVars="0.1 0.05 0.01"
inputDir="data"
execPath="bin/x64/main"
outDir="run/nn_hyperparam_tuning_data"

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
    "33000-86000"
)

declare -a tlim=(
    1
    2
    3.5
    6
    10
    18
    32
    60
    100
    200
)


solverFixedArgs="-m nn --graspType almostbest --loglvl log --round"
outFname="nn_almostbest"

for i in ${!subDirs[@]}; do

    solverArgs="$solverFixedArgs -t ${tlim[$i]}"
    outTemplateFname="$outDir/$outFname"_"${subDirs[$i]}_${tlim[$i]}sec"
    outCostsFname="$outTemplateFname-costs.csv"
    outRuntimesFname="$outTemplateFname-runtimes.csv"
    outIterCountFname="$outTemplateFname-iterCount.csv"
    fullInputPath="$inputDir/${subDirs[$i]}"
    
    echo $BenchmarkPy --execPath $execPath -n $nRuns --solverExtraArgs \"$solverArgs\" --inputDir $fullInputPath --saveCosts $outCostsFname --param2Tune $param2Tune --tuningVars $tuningVars
    python $BenchmarkPy --execPath $execPath -n $nRuns --solverExtraArgs "$solverArgs" --inputDir $fullInputPath --saveCosts $outCostsFname --param2Tune $param2Tune --tuningVars $tuningVars
done

solverFixedArgs="-m nn --graspType random --loglvl log --round"
outFname="nn_random"

for i in ${!subDirs[@]}; do

    solverArgs="$solverFixedArgs -t ${tlim[$i]}"
    outTemplateFname="$outDir/$outFname"_"${subDirs[$i]}_${tlim[$i]}sec"
    outCostsFname="$outTemplateFname-costs.csv"
    outRuntimesFname="$outTemplateFname-runtimes.csv"
    outIterCountFname="$outTemplateFname-iterCount.csv"
    fullInputPath="$inputDir/${subDirs[$i]}"
    
    echo $BenchmarkPy --execPath $execPath -n $nRuns --solverExtraArgs \"$solverArgs\" --inputDir $fullInputPath --saveCosts $outCostsFname --param2Tune $param2Tune --tuningVars $tuningVars
    python $BenchmarkPy --execPath $execPath -n $nRuns --solverExtraArgs "$solverArgs" --inputDir $fullInputPath --saveCosts $outCostsFname --param2Tune $param2Tune --tuningVars $tuningVars
done

solverFixedArgs="-m nn --nnTryall --loglvl log --round"
outFname="nn_tryall_deterministic"

for i in ${!subDirs[@]}; do

    solverArgs="$solverFixedArgs -t ${tlim[$i]}"
    outTemplateFname="$outDir/$outFname"_"${subDirs[$i]}_${tlim[$i]}sec"
    outCostsFname="$outTemplateFname-costs.csv"
    outRuntimesFname="$outTemplateFname-runtimes.csv"
    outIterCountFname="$outTemplateFname-iterCount.csv"
    fullInputPath="$inputDir/${subDirs[$i]}"
    
    echo $BenchmarkPy --execPath $execPath -n $nRuns --solverExtraArgs \"$solverArgs\" --inputDir $fullInputPath --saveCosts $outCostsFname
    python $BenchmarkPy --execPath $execPath -n $nRuns --solverExtraArgs "$solverArgs" --inputDir $fullInputPath --saveCosts $outCostsFname
done
