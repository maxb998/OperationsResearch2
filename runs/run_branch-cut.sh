#! /bin/bash

execPath="bin/exec/main"
inputDir="data"

nRuns=5

outDir="runs/branch-cut"
outFname="branch-cut"

declare -a subDirs=("0-80" "100-200" "220-320" "400-500" "500-800")

declare -a tlims=( 5 30 120 300 400 )


solverFixedArgs="-m branch-cut --metaInit nn --cplexInit vns --graspType almostbest --round --cplexDisableSolPosting --cplexDisableUsercuts"
for i in ${!subDirs[@]}; do
    for file in $inputDir/${subDirs[$i]}/*; do
        for repeatCounter in $(seq 1 $nRuns); do
            seed=$SRANDOM
            outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_base\_$seed.log
            while [ -f $outFullFilename ]; do
                echo "File $outFullFilename already exists. Selecting another seed value..."
                seed=$SRANDOM
                outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_base\_$seed.log
            done
            echo $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed \> $outFullFilename
            $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed > $outFullFilename
        done
    done
done

solverFixedArgs="-m branch-cut --metaInit nn --cplexInit vns --graspType almostbest --round --cplexEnableWarmStart --cplexDisableSolPosting --cplexDisableUsercuts"
for i in ${!subDirs[@]}; do
    for file in $inputDir/${subDirs[$i]}/*; do
        for repeatCounter in $(seq 1 $nRuns); do
            seed=$SRANDOM
            outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_warmstart5\_$seed.log
            while [ -f $outFullFilename ]; do
                echo "File $outFullFilename already exists. Selecting another seed value..."
                seed=$SRANDOM
                outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_warmstart5\_$seed.log
            done
            echo $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed \> $outFullFilename
            $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed > $outFullFilename
        done
    done
done

solverFixedArgs="-m branch-cut --metaInit nn --cplexInit vns --graspType almostbest --round --cplexEnableWarmStart --cplexDisableSolPosting"
for i in ${!subDirs[@]}; do
    for file in $inputDir/${subDirs[$i]}/*; do
        for repeatCounter in $(seq 1 $nRuns); do
            seed=$SRANDOM
            outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_warmstart5+usercuts\_$seed.log
            while [ -f $outFullFilename ]; do
                echo "File $outFullFilename already exists. Selecting another seed value..."
                seed=$SRANDOM
                outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_warmstart5+usercuts\_$seed.log
            done
            echo $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed \> $outFullFilename
            $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed > $outFullFilename
        done
    done
done

solverFixedArgs="-m branch-cut --metaInit nn --cplexInit vns --graspType almostbest --round --cplexEnableWarmStart --cplexDisableUsercuts"
for i in ${!subDirs[@]}; do
    for file in $inputDir/${subDirs[$i]}/*; do
        for repeatCounter in $(seq 1 $nRuns); do
            seed=$SRANDOM
            outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_warmstart5+posting\_$seed.log
            while [ -f $outFullFilename ]; do
                echo "File $outFullFilename already exists. Selecting another seed value..."
                seed=$SRANDOM
                outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_warmstart5+posting\_$seed.log
            done
            echo $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed \> $outFullFilename
            $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed > $outFullFilename
        done
    done
done

solverFixedArgs="-m branch-cut --metaInit nn --cplexInit vns --graspType almostbest --round --cplexEnableWarmStart"
for i in ${!subDirs[@]}; do
    for file in $inputDir/${subDirs[$i]}/*; do
        for repeatCounter in $(seq 1 $nRuns); do
            seed=$SRANDOM
            outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_full\_$seed.log
            while [ -f $outFullFilename ]; do
                echo "File $outFullFilename already exists. Selecting another seed value..."
                seed=$SRANDOM
                outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_full\_$seed.log
            done
            echo $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed \> $outFullFilename
            $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed > $outFullFilename
        done
    done
done