#! /bin/bash

inputDir="data/generated"
execPath="bin/exec/main"

nRuns=10

solverFixedArgs="-m nn --2opt --round"
outDir="runs/nn2opt_MT"
outFname="nn2opt"

param2Tune="threads"
declare -a tuningVars=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 )

declare -a subDirs=("100" "500" "1000" "5000" "10000")
declare -a tlims=( 1 2 4 8 16 )

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