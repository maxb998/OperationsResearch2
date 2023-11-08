import os
import argparse
import subprocess
import numpy as np
from colorama import Fore
from pathlib import Path
import csv
import re
import time 

np.set_printoptions(precision=2)

'''
datadict = {
    'instances', 
    'execPath', 
    'solverExtraArgs', 
    'out_fnames',
    'nIters',
    'param2Tune',
    'tuningVars'
}
'''

def main():

    datadict = arg_parser()

    os.system('clear')

    if len(datadict['param2Tune']) > 0:
        print('Hyperparameter to tune is ' + datadict['param2Tune'] + ' to tune with values: ' + Fore.LIGHTYELLOW_EX + str(datadict['tuningVars']) + Fore.RESET)
    
    print('\n')

    runtimes_table = np.zeros([len(datadict['instances']), len(datadict['tuningVars'])], dtype=float)
    costs_table = np.zeros(runtimes_table.shape, dtype=float)

    for instIndex in range(len(datadict['instances'])):

        print('Running on instance: ' + Fore.WHITE + Path(datadict['instances'][instIndex]).stem + Fore.RESET)

        for tuneValIndex in range(len(datadict['tuningVars'])):

            if len(datadict['param2Tune']) > 0:
                print('\tRunning with ' + datadict['param2Tune'] + '  ' + Fore.LIGHTYELLOW_EX + datadict['tuningVars'][tuneValIndex] + Fore.RESET)

            # could also sum the result for direct average but might need this in the future
            costResults = np.zeros(shape=datadict['nIters'], dtype=float)
            runtimeResults = np.zeros(shape=datadict['nIters'], dtype=float)

            for i in range(datadict['nIters']):

                seed = np.random.randint(low=0, high=2147483648)

                # run the solver on the first instance
                cmd = get_cmd_list(datadict, datadict['instances'][instIndex], datadict['tuningVars'][tuneValIndex], seed)
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

                    elif 'Total runtime = ' in line:
                        start_pos = line.find('=') + 2
                        end_pos = line.find(' seconds',  start_pos)
                        runtimeResults[i] = float(line[ start_pos : end_pos])

                    elif 'Final cost = ' in line:
                        start_pos = line.find('=') + 2
                        end_pos = line.find('\\',  start_pos)
                        costResults[i] = float(line[ start_pos : end_pos])
                
                
                costLogStr = 'Cost = ' + str('{0:.2f}').format(costResults[i])
                runtimeLogStr = 'Runtime = ' + str('{0:.2f}').format(runtimeResults[i]) + ' s'
                seedLogStr = 'Seed = ' + str(seed)
                print('\t\t' + f'{costLogStr: <35}{runtimeLogStr: <30}{seedLogStr}')
            
            runtimes_table[instIndex, tuneValIndex] = np.average(runtimeResults)
            costs_table[instIndex, tuneValIndex] = np.average(costResults)

            costLogStr = 'Avg Cost = ' + Fore.LIGHTGREEN_EX + str('{0:.2f}').format(costs_table[instIndex, tuneValIndex]) + Fore.RESET
            runtimeLogStr = 'Avg Runtime = ' + Fore.LIGHTCYAN_EX + str('{0:.2f}').format(runtimes_table[instIndex, tuneValIndex]) + Fore.RESET + ' s'
            print('\t' + f'{costLogStr: <35}{runtimeLogStr}')
        
    write_csv(runtimes_table, datadict['out_fnames'][0], datadict)
    write_csv(costs_table, datadict['out_fnames'][1], datadict)



def arg_parser() -> dict:

    # parse arguments
    parser = argparse.ArgumentParser(prog='TspBenchmark', epilog='Tsp Solver Benchmarking Tool', description='Run tsp solver multiple times saving execution times, number of iterations and cost of the solutions found.')
    parser.add_argument('--inputDir', metavar='str', required=True, type=str, help='Directory containing all the tsp file instances on which run the solver.')
    parser.add_argument('--execPath', metavar='str', required=True, type=str, help='Location of the exec.')
    parser.add_argument('--outputFile', metavar='str', required=True, type=str, help='Name of the output file')
    parser.add_argument('-n', '--nIters', metavar='int', required=True, type=int, help='Number of times that the runs are repeated with a different random seed')
    parser.add_argument('--param2Tune', metavar='str', type=str, help='Specify an hyperparameter to tune.')
    parser.add_argument('--tuningVars', metavar='str', type=str, nargs='+', help='Specify the values to assign to the hyperparameter specified by param2Tune to use.')
    parser.add_argument('--solverExtraArgs', metavar='str', required=True, type=str, help='Extra commandline arguments to pass to the solver.')

    args = parser.parse_args()
    
    datadict = {}

    # PATHS
    if not os.path.isdir(args.inputDir):
        print('inputDir must be a directory')
        exit()
    datadict['instances'] = list(filter( lambda x: os.path.splitext(x)[1] == '.tsp', os.listdir(args.inputDir) ))

    # sort by number of nodes specified as the first number in the filename
    instSizes = np.zeros(shape=len(datadict["instances"]), dtype=np.int32)
    for i in range(len(datadict["instances"])):
        instSizes[i] = re.findall(r'\d+', datadict["instances"][i])[0]
    sortedIndexes = np.argsort(instSizes)
    datadict["instances"] = [datadict["instances"][i] for i in sortedIndexes]

    datadict['instances'] = list( filter( lambda x: os.path.isfile(os.path.join(args.inputDir, x)), iter(datadict['instances']) ))
    for i in range(len(datadict['instances'])):
        datadict['instances'][i] = os.path.join(args.inputDir, datadict['instances'][i])

    if (len(datadict['instances']) == 0):
        print('There are no compatible files in inputDir')
        exit()

    if args.execPath == None or not os.path.isfile(args.execPath):
        print('execPath is not detected/not a file')
        exit()
    datadict['execPath'] = args.execPath

    datadict['out_fnames'] = [args.outputFile + '_costs.csv', args.outputFile + 'runtimes.csv']
    for out_fname in datadict['out_fnames']:
        if os.path.isfile(out_fname):
            while True:
                print('File ' + out_fname + ' already exists. Overwrite? (y/n)')
                c = input('')[0]
                print ('\033[A                             \033[A')
                if c == 'y':
                    print('File will be overwritten')
                    break
                elif c == 'n':
                    print('File will not be overwritten. Quitting...')
                    exit()
                print ('\033[A                             \033[A')

    datadict['nIters'] = args.nIters

    # get args for tsp solver
    datadict['solverExtraArgs'] = args.solverExtraArgs.split(' ')
    
    # get hyperparameter tuning values
    if args.param2Tune != None:
        if args.tuningVars == None:
            print(Fore.LIGHTRED_EX + 'Must specify --tuningVars along with --param2Tune with the desired values to test' + Fore.RESET)
            exit()
        datadict['param2Tune'] = args.param2Tune
        datadict['tuningVars'] = np.array(args.tuningVars)
    else:
        datadict['param2Tune'] = ['']
        datadict['tuningVars'] = ['']

    return datadict


def get_cmd_list(datadict:dict, inst:str, val:str, seed:int) -> list[str]:
    cmd = []
    cmd.append(datadict['execPath'])
    for solver_arg in datadict['solverExtraArgs']:
        cmd.append(solver_arg)
    cmd.append('-f')
    cmd.append(inst)
    cmd.append('--seed')
    cmd.append(str(seed))
    cmd.append('--loglvl')
    cmd.append('notice')
    if len(datadict['param2Tune']) > 0:
        cmd.append('--' + datadict['param2Tune'])
        cmd.append(val)
    return cmd


def write_csv(table:np.ndarray, csv_filename:str, datadict:dict):

    with open(csv_filename, mode='a', newline='') as csvfile:

        csv_writer = csv.writer(csvfile, delimiter=' ')

        lineList = [ str(len(datadict['tuningVars'])) ]
        for paramVal in datadict['tuningVars']:
            lineList.append(datadict['param2Tune'] + '=' + paramVal)

        csv_writer.writerow(lineList)

        for instance_num in range(table.shape[0]):
            lineList.clear()
            lineList.append(Path(datadict['instances'][instance_num]).stem)
            for iter_num in range(table.shape[1]):
                lineList.append(str(table[instance_num, iter_num]))
            csv_writer.writerow(lineList)



if __name__ == '__main__' :
    main()
