#! /bin/bash

execPath="bin/x64/main"
inputDir="data"

nRuns=3

solverFixedArgs="-m tabu --metaInit nn --graspType almostbest --round"
outDir="runs/tabu"
outFname="tabu"

param2Tune="graspChance"
declare -a tabuTenureSize=( 5 20 60 )

declare -a subDirs=("0-80" "100-200" "220-320" "400-500" "500-800" "1000-1440" "1570-2400" "3000-6000")

declare -a tlims=  ( 1      3         6         20        60        180         400         600       )

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

