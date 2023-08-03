import os
import argparse
import subprocess
import random
from datetime import datetime
import numpy as np
import csv

'''
datadict = {
    "instances", 
    "exec", 
    "solver_args", 
    "out_fnames", 
    "num_iters",
    "mode",
    "seed"
}
'''

#from matplotlib import pyplot

#seeds_list = [ 2023, 567441, 25, 789632175, 5, 855699, 111111]

def main():
    random.seed(datetime.now().timestamp())

    datadict = arg_parser()

    datadict["seed"] = np.zeros(datadict["num_iters"], dtype=np.int32)

    runtimes_table = np.zeros([datadict["num_iters"], len(datadict["instances"])], dtype=float)
    costs_table = np.zeros(runtimes_table.shape, dtype=float)

    for i in range(datadict["num_iters"]):

        datadict["seed"][i] = random.randint(0, 2147483647)
        print("Iteration seed = " + str(datadict["seed"][i]))

        for j in range(len(datadict["instances"])):

            # run the solver on the first instance
            cmd = get_cmd_list(datadict, j, i)
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

            exec_out = []

            for line in p.stdout:
                line = str(line)
                exec_out.append(line)

                if "ERR" in line:
                    for l in exec_out:
                        print(l)
                    exit(1)
                
                elif "WARN" in line:
                    print(line)

                elif "Total runtime = " in line:
                    start_pos = line.find("=") + 2
                    end_pos = line.find(" seconds",  start_pos)
                    runtimes_table[i,j] = float(line[ start_pos : end_pos])

                elif "Final cost = " in line:
                    start_pos = line.find("=") + 2
                    end_pos = line.find("\\",  start_pos)
                    costs_table[i,j] = float(line[ start_pos : end_pos])
            
            print("   Run for " + str(datadict["instances"][j]) + " has finished with cost " + str(costs_table[i,j]) + "   in " + str(runtimes_table[i,j]) + "seconds")
    
    write_csv(runtimes_table, datadict["out_fnames"][0], datadict)
    write_csv(costs_table, datadict["out_fnames"][1], datadict)



def arg_parser() -> dict:

    # parse arguments
    parser = argparse.ArgumentParser(prog='TspBenchmark', epilog='Tsp Solver Benchmarking Tool', description='Run tsp solver multiple times saving execution times, number of iterations and cost of the solutions found')
    parser.add_argument('input_dir', type=str, help='Directory containing all the tsp file instances on which run the solver')
    parser.add_argument('exec', type=str, help='Location of the exec')
    parser.add_argument('solver_args', type=str, help='Commandline arguments to pass to the solver')
    parser.add_argument('output_filename', type=str, help='Name of the output file')
    parser.add_argument('-n', '--num_iters', metavar='ITERNUM', type=int, help='Number of runs per-instance')

    args = parser.parse_args()
    
    datadict = {}

    if not os.path.isdir(args.input_dir):
        print('input_dir must be a directory')
        exit()
    datadict["instances"] = filter( lambda x: os.path.splitext(x)[1] == '.tsp', os.listdir(args.input_dir) )
    datadict["instances"] = sorted( filter( lambda x: os.path.isfile(os.path.join(args.input_dir, x)), iter(datadict["instances"]) ))
    if (len(datadict["instances"]) == 0):
        print('There are no compatible files in input_dir')
        exit()

    if not os.path.isfile(args.exec):
        print('exec is not detected/not a file')
        exit()
    datadict["exec"] = args.exec

    # check solver_args? no need but would be nice
    datadict["solver_args"] = args.solver_args.split(' ')

    for i in range(len(datadict["solver_args"])):
        if "-m" in datadict["solver_args"][i] or "--mode" in datadict["solver_args"][i]:
            datadict["mode"] = datadict["solver_args"][i+1]

    datadict["out_fnames"] = [args.output_filename + "_costs.csv", args.output_filename + "runtimes.csv"]
    for out_fname in datadict["out_fnames"]:
        if os.path.isfile(out_fname):
            while True:
                print('File ' + out_fname + ' already exists. Overwrite? (y/n)')
                c = input("")[0]
                print ("\033[A                             \033[A")
                if c == 'y':
                    print('File will be overwritten')
                    break
                elif c == 'n':
                    print('File will not be overwritten. Quitting...')
                    exit()
                print ("\033[A                             \033[A")

    datadict["num_iters"]: int = 1
    if args.num_iters != None:
        '''if args.num_iters <= 0 or args.num_iters > len(seeds_list):
            print("Invalid number of iterations value. Must be between range(0, " + str(len(seeds_list)) + ")")
            exit()'''
        datadict["num_iters"]: int = args.num_iters
    
    return datadict


def get_cmd_list(datadict:dict, instance_index:int, iteration:int) -> list[str]:
    cmd = []
    cmd.append(datadict["exec"])
    for solver_arg in datadict["solver_args"]:
        cmd.append(solver_arg)
    cmd.append("-f")
    cmd.append(datadict["instances"][instance_index])
    cmd.append("--seed")
    cmd.append(str(datadict["seed"][iteration]))
    cmd.append("--loglvl")
    cmd.append("notice")
    return cmd


def write_csv(table:np.ndarray, csv_filename:str, datadict:dict):

    with open(csv_filename, mode='w', newline='') as csvfile:

        csv_writer = csv.writer(csvfile, delimiter=' ')

        csv_writer.writerow([ str(1) , datadict["mode"] ])

        for instance_num in range(table.shape[1]):
            for iter_num in range(table.shape[0]):
                csv_writer.writerow([ datadict["instances"][instance_num] + "_SEED=" + str(datadict["seed"][iter_num]) , str(table[iter_num, instance_num]) ])



if __name__ == "__main__" :
    main()
