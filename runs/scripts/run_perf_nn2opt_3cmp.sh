#! /bin/bash

inputDir="data/generated"

nRuns=10

solverFixedArgs="-m nn --round --2opt -j 1"
outDir="runs/nn2opt_3cmp"
outFname="nn2opt"

declare -a execPaths=("bin/exec/main_avx_nocc_noapprox" "bin/exec/main_avx_nocc" "bin/exec/main_avx_noapprox" "bin/exec/main_avx" "bin/exec/main_base_nocc" "bin/exec/main_base" "bin/exec/main_matrix_nocc" "bin/exec/main_matrix_nocc_badindex" "bin/exec/main_matrix")
declare -a mode=("avx-nocc-noapprox" "avx-nocc" "avx-noapprox" "avx" "base-nocc" "base" "matrix-nocc" "matrix-nocc-badindex" "matrix")

declare -a subDirs=("100" "500" "1000" "5000" "10000")
declare -a tlims=( 1 2 4 8 16 )

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
