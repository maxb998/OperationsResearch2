#! /bin/bash

execPath="bin/exec/main"
inputDir="data"

nRuns=3

solverFixedArgs="-m genetic --metaInit nn --graspType almostbest --round"
outDir="runs/genetic"
outFname="genetic"

param2Tune="geneticParams"
declare -a tuningVars=( "100,50,20,20" "100,20,50,20" "100,40,40,10" "50,20,20,5" "50,10,10,5", "100,100,100,100" )

declare -a subDirs=("0-80" "100-200" "220-320" "400-500")

declare -a tlims=( 20 30 60 100 )

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

