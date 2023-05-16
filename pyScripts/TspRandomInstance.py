import os
import argparse
import numpy as np
from matplotlib import pyplot
from datetime import datetime

def main():
    n, dist_type, filename, seed = arg_parser()

    f = open(filename, "w")
    
    f.write("NAME : {0}\n".format(os.path.splitext(os.path.basename(filename))[0]))
    f.write("COMMENT : Randomly generated tsp instance using seed = {0}\n".format(seed))
    f.write("TYPE : TSP\n")
    f.write("DIMENSION : {0}\n".format(n))
    f.write("EDGE_WEIGHT_TYPE : {0}\n".format(dist_type))
    f.write("NODE_COORD_SECTION\n")

    # set seed
    np.random.seed(seed)

    # random points within range [0, 1)
    points = np.random.rand(n, 2)

    # set the points to have range [0, 10000)
    points = points * 10000

    for i in range(n):
        f.write("{0} {1} {2}\n".format(i+1, points[i,0], points[i,1]))

    f.write("EOF\n")

    f.close()

    # plot istance
    x, y = np.hsplit(points, 2)
    pyplot.scatter(x, y, color='red')
    pyplot.title("Instance")
    pyplot.show()
    

def arg_parser() -> tuple[ int, str, str, int]:
    # parse arguments
    parser = argparse.ArgumentParser(prog='TspRandomInstance', epilog='Tsp Random Instance Generator', description='Generate random instance of a traveling salesman problem and saves it according to tsplib specifications.')
    parser.add_argument('number_of_nodes', type=int, help='Number of nodes that the newly generated instance will contain.')
    parser.add_argument('distance_type', choices=["EUC_2D","MAN_2D","MAX_2D", "ATT"], help='Type of distance used. Supported types only.')
    parser.add_argument('-o', '--output_file', metavar='OUTPUT_FILENAME')
    parser.add_argument('-s', '--seed', metavar='SEED', type=int, help='Set random seed manually. (use time in seconds otherwise)')

    args = parser.parse_args()
    
    n: int = args.number_of_nodes
    if n <= 0:
        print('Invalid number of nodes')
        exit()

    dist_type = args.distance_type

    # get random seed
    seed = int(round(datetime.now().timestamp()))
    if args.seed != None:
        if args.seed <= 0:
            print('Invalid seed value')
            exit()
        seed: int = args.seed
    
    # set output file
    if args.output_file == None:
        filename = "data/{0}-{1}-{2}.tsp".format(n, dist_type, seed)
    else:
        filename = args.output_file
    
    if os.path.isfile(filename):
        while True:
            print('File ' + filename + ' already exists. Overwrite? (y/n)')
            c = input("")[0]
            print ("\033[A                             \033[A")
            if c == 'y':
                print('Overwriting...')
                break
            elif c == 'n':
                print('File will not be overwritten. Quitting...')
                exit()
            print ("\033[A                             \033[A")

    return [n, dist_type, filename, seed]



if __name__ == "__main__" :
    main()