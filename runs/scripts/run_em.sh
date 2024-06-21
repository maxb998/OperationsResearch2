#! /bin/bash

execPath="bin/exec/main"
inputDir="data"

solverFixedArgs="-m em --round"
outDir="runs/em/em_nograsp"
outFname="em_nograsp"

nRuns=3

# NOGRASP #################################################################################################################################################

declare -a subDirs=("0-80" "100-200" "220-320" "400-500" "500-800" "1000-1440" "1570-2400")
declare -a tlims=( 1 2 3.5 6 10 18 32 )

for i in ${!subDirs[@]}; do
    for file in $inputDir/${subDirs[$i]}/*; do
        seed=$SRANDOM
        outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_0\_$seed.log
        while [ -f $outFullFilename ]; do
            echo "File $outFullFilename already exists. Selecting another seed value..."
            seed=$SRANDOM
            outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_0\_$seed.log
        done
        echo $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed \> $outFullFilename
        $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed > $outFullFilename
    done
done

# ALMOSTBEST #################################################################################################################################################

solverFixedArgs="-m em --graspType almostbest --round"
outDir="runs/em/em_almostbest"
outFname="em_almostbest"
param2Tune="graspChance"

declare -a tuningVars=( "0.5", "0.4", "0.3", "0.2", "0.15", "0.1", "0.075", "0.05", "0.03", "0.01", "0.001" )
declare -a subDirs=("0-80" "100-200" "220-320" "400-500" "500-800" "1000-1440" "1570-2400" "3000-6000" "7000-20000" "33000-86000")
declare -a tlims=( 1 2 3.5 6 10 18 32 60 100 200 )

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

# RANDOM #################################################################################################################################################

solverFixedArgs="-m em --graspType random --round"
outDir="runs/em/em_random"
outFname="em_random"
param2Tune="graspChance"

declare -a tuningVars=( "0.1", "0.05", "0.03", "0.01", "0.001" )
declare -a subDirs=("0-80" "100-200" "220-320" "400-500" "500-800" "1000-1440" "1570-2400")
declare -a tlims=( 1 2 3.5 6 10 18 32 )

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

declare -a tuningVars=( "0.02","0.005" )
declare -a subDirs=("0-80" "100-200" "220-320" "400-500" "500-800")
declare -a tlims=( 1 2 3.5 6 10 )

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