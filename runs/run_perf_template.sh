#! /bin/bash

inputDir="data/generated"

# nRuns=3

# solverFixedArgs="-m nn --round --2opt -j 1"
# outDir="runs/more_avxshowcase"
# outFname="nn2opt"

# declare -a execPaths=("bin/exec/main_avx" "bin/exec/main_base" "bin/exec/main_matrix")
# declare -a mode=("avx" "base" "matrix")

# declare -a subDirs=("100" "500" "1000" "5000" "10000")
# declare -a tlims=( 3 10 30 90 270)

# for i in ${!subDirs[@]}; do
#     for file in $inputDir/${subDirs[$i]}/*; do
#         for j in ${!execPaths[@]}; do
#             for repeatCounter in $(seq 1 $nRuns); do
#                 seed=$SRANDOM
#                 outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_${mode[$j]}\_$seed.log
#                 while [ -f $outFullFilename ]; do
#                     echo "File $outFullFilename already exists. Selecting another seed value..."
#                     seed=$SRANDOM
#                     outFullFilename=$outDir/$outFname\_$(basename ${file%.*})\_${tlims[$i]}s\_${mode[$j]}\_$seed.log
#                 done
#                 echo ${execPaths[$j]} -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed \> $outFullFilename
#                 ${execPaths[$j]} -f $file -t ${tlims[$i]} $solverFixedArgs --seed $seed > $outFullFilename
#             done
#         done
#     done
# done

execPath="bin/exec/main"

nRuns=1

solverFixedArgs="-m nn --2opt --round"
outDir="runs/more_MTshowcase"
outFname="nn2opt"

param2Tune="threads"
declare -a tuningVars=( 1 )

declare -a subDirs=("100" "500" "1000" "5000" "10000")
declare -a tlims=( 3 10 30 90 270)

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