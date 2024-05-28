#! /bin/bash

inputDir="data/perf"

nRuns=3

solverFixedArgs="-m nn --round --2opt"
outDir="runs/perf"
outFname="nn2opt"

declare -a execPaths=("bin/exec/main_avx" "bin/exec/main_base" "bin/exec/main_matrix")
declare -a mode=("AVX" "BASE" "MATRIX")

declare -a subDirs=("100" "500" "1000" "5000")
declare -a tlims=( 3 10 30 90 )

for i in ${!subDirs[@]}; do
    for file in $inputDir/${subDirs[$i]}/*; do
        for j in ${!execPaths[@]}; do
            for repeatCounter in $(seq 1 $nRuns); do
                seed=$SRANDOM
                outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_${mode[$j]}\_$seed.log
                while [ -f $outFullFilename ]; do
                    echo "File $outFullFilename already exists. Selecting another seed value..."
                    seed=$SRANDOM
                    outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_${mode[$j]}\_$seed.log
                done
                echo ${execPaths[$j]} -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed \> $outFullFilename
                ${execPaths[$j]} -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed > $outFullFilename
            done
        done
    done
done

