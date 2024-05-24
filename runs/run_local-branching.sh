#! /bin/bash

execPath="bin/exec/main"
inputDir="data"

nRuns=5

outDir="runs/local-branching"
outFname="local-branching"

declare -a subDirs=("0-80" "100-200" "220-320" "400-500" "500-800" "1000-1440" "1570-2400")

declare -a tlims=( 1 3 8 20 60 180 400)

solverFixedArgs="-m local-branching --metaInit nn --cplexInit vns --graspType almostbest --round "
for i in ${!subDirs[@]}; do
    for file in $inputDir/${subDirs[$i]}/*; do
        for repeatCounter in $(seq 1 $nRuns); do
            seed=$SRANDOM
            outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_10-5\_$seed.log
            while [ -f $outFullFilename ]; do
                echo "File $outFullFilename already exists. Selecting another seed value..."
                seed=$SRANDOM
                outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_10-5\_$seed.log
            done
            echo $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed \> $outFullFilename
            $execPath -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed > $outFullFilename
        done
    done
done
