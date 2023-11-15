#! /bin/bash

BenchmarkPy="pyScripts/Benchmark.py"

nRuns=3
param2Tune="graspChance"
tuningVars="0.1 0.05 0.01"
inputDir="data"
execPath="bin/x64/main"
solverFixedArgs="-m em --graspType almostbest --loglvl log --round"
outDir="run"
outFname="em_almostbest"

declare -a subDirs=(
    "0-80"
    "100-200"
    "220-320"
    "400-500"
    "500-800"
    "1000-1440"
    "1570-2400"
    "3000-6000"
    #"7000-20000"
    #"33000-86000"
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
    #100
    #200
)

# should take 5 hours on a 12700k

for i in ${!subDirs[@]}; do

    solverArgs="$solverFixedArgs -t ${tlim[$i]}"
    outTemplateFname="$outDir/$outFname"_"${subDirs[$i]}_${tlim[$i]}sec"
    outCostsFname="$outTemplateFname-costs.csv"
    outRuntimesFname="$outTemplateFname-runtimes.csv"
    outIterCountFname="$outTemplateFname-iterCount.csv"
    fullInputPath="$inputDir/${subDirs[$i]}"
    
    echo $BenchmarkPy --execPath $execPath -n $nRuns --solverExtraArgs \"$solverArgs\" --inputDir $fullInputPath --saveCosts $outCostsFname --saveIterCount $outIterCountFname --param2Tune $param2Tune --tuningVars $tuningVars # --saveRuntimes $outRuntimesFname
    python $BenchmarkPy --execPath $execPath -n $nRuns --solverExtraArgs "$solverArgs" --inputDir $fullInputPath --saveCosts $outCostsFname --saveIterCount $outIterCountFname --param2Tune $param2Tune --tuningVars $tuningVars # --saveRuntimes $outRuntimesFname
done

solverFixedArgs="-m em --graspType random --loglvl log --round"
outFname="em_random"

for i in ${!subDirs[@]}; do

    solverArgs="$solverFixedArgs -t ${tlim[$i]}"
    outTemplateFname="$outDir/$outFname"_"${subDirs[$i]}_${tlim[$i]}sec"
    outCostsFname="$outTemplateFname-costs.csv"
    outRuntimesFname="$outTemplateFname-runtimes.csv"
    outIterCountFname="$outTemplateFname-iterCount.csv"
    fullInputPath="$inputDir/${subDirs[$i]}"
    
    echo $BenchmarkPy --execPath $execPath -n $nRuns --solverExtraArgs \"$solverArgs\" --inputDir $fullInputPath --saveCosts $outCostsFname --saveIterCount $outIterCountFname --param2Tune $param2Tune --tuningVars $tuningVars # --saveRuntimes $outRuntimesFname
    python $BenchmarkPy --execPath $execPath -n $nRuns --solverExtraArgs "$solverArgs" --inputDir $fullInputPath --saveCosts $outCostsFname --saveIterCount $outIterCountFname --param2Tune $param2Tune --tuningVars $tuningVars # --saveRuntimes $outRuntimesFname
done

solverFixedArgs="-m em --emFarthest --loglvl log --round"
outFname="em_fathest_deterministic"

for i in ${!subDirs[@]}; do

    solverArgs="$solverFixedArgs -t ${tlim[$i]}"
    outTemplateFname="$outDir/$outFname"_"${subDirs[$i]}_${tlim[$i]}sec"
    outCostsFname="$outTemplateFname-costs.csv"
    outRuntimesFname="$outTemplateFname-runtimes.csv"
    outIterCountFname="$outTemplateFname-iterCount.csv"
    fullInputPath="$inputDir/${subDirs[$i]}"
    
    echo $BenchmarkPy --execPath $execPath -n $nRuns --solverExtraArgs \"$solverArgs\" --inputDir $fullInputPath --saveCosts $outCostsFname --saveIterCount $outIterCountFname # --saveRuntimes $outRuntimesFname
    python $BenchmarkPy --execPath $execPath -n $nRuns --solverExtraArgs "$solverArgs" --inputDir $fullInputPath --saveCosts $outCostsFname --saveIterCount $outIterCountFname # --saveRuntimes $outRuntimesFname
done

#python pyScripts/Benchmark.py -n 1 --solverExtraArgs "-m nn -t 2" --inputDir data/0-80/ --saveCosts run/testCost.csv --saveIterCount run/testIterCount.csv --saveRuntimes run/testRuntimes.csv