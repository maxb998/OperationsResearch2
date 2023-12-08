#! /bin/bash

execPath="bin/x64/main"
inputDir="data"

nRuns=8

solverFixedArgs="-m em --graspType almostbest --loglvl log --round"
outDir="run/em_almostbest"
outFname="em_almostbest"

param2Tune="graspChance"
declare -a tuningVars=(
    0.10
    0.05
    0.01
)

declare -a subDirs=(
    "0-80"
    "100-200"
    "220-320"
    "400-500"
    "500-800"
    "1000-1440"
    "1570-2400"
    "3000-6000"
    # "7000-20000"
    # "33000-86000"
)

declare -a tlims=(
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

for i in ${!subDirs[@]}; do
    for file in $inputDir/${subDirs[$i]}/*; do
        for j in ${!tuningVars[@]}; do
            for repeatCounter in $(seq 1 $nRuns); do
                seed=$SRANDOM
                outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_${tuningVars[$j]}\_$seed.log
                while [ -f $outFullFilename ]; do
                    echo "File $outFullFilename already exists. Selecting another seed value..."
                    seed=$SRANDOM
                    outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_${tuningVars[$j]}\_$seed.log
                done
                echo $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --$param2Tune ${tuningVars[$j]} --seed $seed \> $outFullFilename
                $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --$param2Tune ${tuningVars[$j]} --seed $seed > $outFullFilename
            done
        done
    done
done


solverFixedArgs="-m em --graspType random --loglvl log --round"
outDir="run/em_random"
outFname="em_random"


for i in ${!subDirs[@]}; do
    for file in $inputDir/${subDirs[$i]}/*; do
        for j in ${!tuningVars[@]}; do
            for repeatCounter in $(seq 1 $nRuns); do
                outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_${tuningVars[$j]}\_$SRANDOM.log
                while [ -f $outFullFilename ]; do
                    outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_${tuningVars[$j]}\_$SRANDOM.log
                done
                echo $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --$param2Tune ${tuningVars[$j]} --seed $seed \> $outFullFilename
                $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --$param2Tune ${tuningVars[$j]} --seed $seed > $outFullFilename
            done
        done
    done
done
