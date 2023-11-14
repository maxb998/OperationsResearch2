import os
import argparse
import subprocess
import numpy as np
from colorama import Fore
from pathlib import Path
import csv
import re
import time 
import random

np.set_printoptions(precision=2)


class Params:
    instances = []
    execPath = 'bin/x64/main'
    solverExtraArgs = ''
    nIters = 10
    param2Tune = ''
    tuningVars = ['']
    costFilename = ''
    runtimesFilename = ''
    iterCountFilename = ''
    

def main():

    random.seed(time.time())

    params = arg_parser()

    if params.param2Tune != '':
        print('Hyperparameter to tune is ' + params.param2Tune + ' to tune with values: ' + Fore.LIGHTYELLOW_EX + str(params.tuningVars) + Fore.RESET)
    
    print("======================================================================================================================================================")

    costs_table = np.zeros([len(params.instances), len(params.tuningVars)], dtype=float)
    runtimes_table = np.zeros(costs_table.shape, dtype=float)
    iterCount_table = np.zeros(costs_table.shape, dtype=float)

    for instIndex in range(len(params.instances)):

        for tuneValIndex in range(len(params.tuningVars)):

            print('Running on instance: ' + Fore.WHITE + Path(params.instances[instIndex]).stem + Fore.RESET + ' ...')

            # could also sum the result for direct average but might need this in the future
            costResults = np.zeros(shape=params.nIters, dtype=float)
            runtimeResults = np.zeros(shape=params.nIters, dtype=float)
            iterCountResult = np.zeros(shape=params.nIters, dtype=float)

            for i in range(params.nIters):

                seed = random.randint(0, 2147483647)

                # run the solver on the first instance
                cmd = get_cmd_list(params, params.instances[instIndex], params.tuningVars[tuneValIndex], seed)
                #print(cmd)
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

                exec_out = []

                for line in p.stdout:
                    line = str(line)
                    exec_out.append(line)

                    if 'ERR' in line:
                        for l in exec_out:
                            print(l)
                        exit(1)
                    
                    elif 'WARN' in line:
                        print(line)

                    elif 'Final cost = ' in line:
                        start_pos = line.find('=') + 2
                        end_pos = line.find('\\',  start_pos)
                        costResults[i] = float(line[ start_pos : end_pos])
                        
                    elif 'finished in ' in line:
                        start_pos = line.find('finished in ') + 12
                        end_pos = line.find(' second',  start_pos)
                        runtimeResults[i] = float(line[ start_pos : end_pos])

                    elif 'Iterations-per-second' in line:
                        start_pos = line.find(':') + 2
                        end_pos = line.find('\\',  start_pos)
                        iterCountResult[i] = float(line[ start_pos : end_pos])
            
            costs_table[instIndex, tuneValIndex] = np.average(costResults)
            runtimes_table[instIndex, tuneValIndex] = np.average(runtimeResults)
            iterCount_table[instIndex, tuneValIndex] = np.average(iterCountResult)

            print ("\033[A                             \033[A") #clear last line of output
            instanceLogStr = 'Result on instance: ' + Fore.WHITE + Path(params.instances[instIndex]).stem + Fore.RESET
            param2TuneLogStr = ''
            if params.param2Tune != '':
                param2TuneLogStr = params.param2Tune + ' = ' + Fore.LIGHTYELLOW_EX + params.tuningVars[tuneValIndex] + Fore.RESET
            costLogStr = 'Avg Cost = ' + Fore.LIGHTGREEN_EX + str('{0:.2f}').format(costs_table[instIndex, tuneValIndex]) + Fore.RESET
            runtimeLogStr = 'Avg Runtime = ' + Fore.LIGHTCYAN_EX + str('{0:.2f}').format(runtimes_table[instIndex, tuneValIndex]) + Fore.RESET + ' s'
            iterCountLogStr = 'Avg Iters/sec = ' + Fore.LIGHTBLUE_EX + str('{0:11.0f}').format(iterCount_table[instIndex, tuneValIndex]) + Fore.RESET + ' iter/s'
            print(f'{instanceLogStr: <51}{param2TuneLogStr: <35}{costLogStr: <35}{runtimeLogStr: <35}{iterCountLogStr}')
            
        print("======================================================================================================================================================")
    
    if params.costFilename != '':
        write_csv(costs_table, params.costFilename, params)
    if params.runtimesFilename != '':
        write_csv(runtimes_table, params.runtimesFilename, params)
    if params.iterCountFilename != '':
        write_csv(iterCount_table, params.iterCountFilename, params)


def arg_parser() -> Params:

    # parse arguments
    parser = argparse.ArgumentParser(prog='TspBenchmark', epilog='Tsp Solver Benchmarking Tool', description='Run tsp solver multiple times saving execution times, number of iterations and cost of the solutions found.')
    parser.add_argument('--inputDir', metavar='str', required=True, type=str, help='Directory containing all the tsp file instances on which run the solver.')
    parser.add_argument('--execPath', metavar='str', required=False, type=str, help='Location of the exec.')
    parser.add_argument('-n', '--nIters', metavar='int', required=True, type=int, help='Number of times that the runs are repeated with a different random seed')
    parser.add_argument('--param2Tune', metavar='str', type=str, help='Specify an hyperparameter to tune.')
    parser.add_argument('--tuningVars', metavar='str', type=str, nargs='+', help='Specify the values to assign to the hyperparameter specified by param2Tune to use.')
    parser.add_argument('--solverExtraArgs', metavar='str', required=True, type=str, help='Extra commandline arguments to pass to the solver.')
    parser.add_argument('--saveCosts', metavar='str', required=False, type=str, help='Output filename for the costs csv')
    parser.add_argument('--saveRuntimes', metavar='str', required=False, type=str, help='Output filename for the runtimes csv')
    parser.add_argument('--saveIterCount', metavar='str', required=False, type=str, help='Output filename for the IterCount csv')

    args = parser.parse_args()
    
    params = Params()

    # PATHS
    if not os.path.isdir(args.inputDir):
        print('inputDir must be a directory')
        exit()
    params.instances = list(filter( lambda x: os.path.splitext(x)[1] == '.tsp', os.listdir(args.inputDir) ))

    # sort by number of nodes specified as the first number in the filename
    instSizes = np.zeros(shape=len(params.instances), dtype=np.int32)
    for i in range(len(params.instances)):
        instSizes[i] = re.findall(r'\d+', params.instances[i])[0]
    sortedIndexes = np.argsort(instSizes)
    params.instances = [params.instances[i] for i in sortedIndexes]

    params.instances = list( filter( lambda x: os.path.isfile(os.path.join(args.inputDir, x)), iter(params.instances) ))
    for i in range(len(params.instances)):
        params.instances[i] = os.path.join(args.inputDir, params.instances[i])

    if (len(params.instances) == 0):
        print('There are no compatible files in inputDir')
        exit()

    if args.execPath != None:
        params.execPath = args.execPath
    if not os.path.isfile(params.execPath):
        print('execPath is not detected/not a file')
        exit()

    if args.saveCosts != None:
        if not os.path.isdir(os.path.dirname(args.saveCosts)):
            print('The directory in which the file for the costs is supposed to be saved does not exists')
            exit()
        params.costFilename = args.saveCosts

    if args.saveRuntimes != None:
        if not os.path.isdir(os.path.dirname(args.saveRuntimes)):
            print('The directory in which the file for the costs is supposed to be saved does not exists')
            exit()
        params.runtimesFilename = args.saveRuntimes
    
    if args.saveIterCount != None:
        if not os.path.isdir(os.path.dirname(args.saveIterCount)):
            print('The directory in which the file for the costs is supposed to be saved does not exists')
            exit()
        params.iterCountFilename = args.saveIterCount


    if (os.path.isfile(params.costFilename) and params.costFilename != '') or (os.path.isfile(params.runtimesFilename) and params.runtimesFilename != '') or (os.path.isfile(params.iterCountFilename) and params.iterCountFilename != ''):
        while True:
            print('One or more output csv file already exists. Overwrite all? (y/n) ', end='')
            c = input('')[0]
            if c == 'y':
                break
            elif c == 'n':
                print('Quitting...')
                exit()

    params.nIters = args.nIters

    # get args for tsp solver
    params.solverExtraArgs = args.solverExtraArgs.split(' ')
    
    # get hyperparameter tuning values
    if args.param2Tune != None:
        if args.tuningVars == None:
            print(Fore.LIGHTRED_EX + 'Must specify --tuningVars along with --param2Tune with the desired values to test' + Fore.RESET)
            exit()
        params.param2Tune = args.param2Tune
        params.tuningVars = np.array(args.tuningVars)

    return params


def get_cmd_list(params:Params, inst:str, val:str, seed:int) -> list[str]:
    cmd = []
    cmd.append(params.execPath)
    for solver_arg in params.solverExtraArgs:
        cmd.append(solver_arg)
    cmd.append('-f')
    cmd.append(inst)
    cmd.append('--seed')
    cmd.append(str(seed))
    cmd.append('--loglvl')
    cmd.append('notice')
    if len(params.param2Tune) > 0:
        cmd.append('--' + params.param2Tune)
        cmd.append(val)
    return cmd


def write_csv(table:np.ndarray, csv_filename:str, params:dict):

    if os.path.isfile(csv_filename):
        os.remove(csv_filename)

    with open(csv_filename, mode='a', newline='') as csvfile:

        csv_writer = csv.writer(csvfile, delimiter=' ')

        lineList = [ str(len(params.tuningVars)) ]
        for paramVal in params.tuningVars:
            lineList.append(params.param2Tune + '=' + paramVal)

        csv_writer.writerow(lineList)

        for instance_num in range(table.shape[0]):
            lineList.clear()
            lineList.append(Path(params.instances[instance_num]).stem)
            for iter_num in range(table.shape[1]):
                lineList.append(str(table[instance_num, iter_num]))
            csv_writer.writerow(lineList)



if __name__ == '__main__' :
    main()
