import sys
import os
import numpy as np
import math
from matplotlib import pyplot


def main():
    hasOptTour = args_check()
    
    x, y, edge_wgt_type = get_points_from_tsp_file(sys.argv[1])

    if edge_wgt_type == -1 :
        print("EDGE_WEIGHT_TYPE of tsp file is not supported")
        exit()

    sol = get_solution_from_tour_file(sys.argv[2])

    if x.shape != sol.shape:
        print("ERROR: tsp file has " + str(x.shape[0]) + " entries while solution file has " + str(sol.shape[0]))
        exit()

    if not hasOptTour:
        x = x[sol[::-1]]
        y = y[sol[::-1]]

        cost = compute_total_cost(x, y, edge_wgt_type)

        pyplot.plot(x, y, color='red', marker='.', linestyle='dashed')
        pyplot.title("Solution with cost: " + str(cost))
        pyplot.show()
    
    else:
        opt_sol = get_solution_from_tour_file(sys.argv[3])

        if x.shape != opt_sol.shape:
            print("ERROR: tsp file has " + str(x.shape[0]) + " entries while the optimal solution file has " + str(opt_sol.shape[0]))
            exit()

        opt_x = x[opt_sol[::-1]]
        opt_y = y[opt_sol[::-1]]

        cost_opt = compute_total_cost(opt_x, opt_y, edge_wgt_type)

        x = x[sol[::-1]]
        y = y[sol[::-1]]

        cost = compute_total_cost(x, y, edge_wgt_type)

        '''f1 = pyplot.figure(1)
        pyplot.plot(x, y, color='red', marker='.', linestyle='dashed')
        pyplot.title("Solution with cost: " + str(cost))

        f2 = pyplot.figure(2)
        pyplot.plot(opt_x, opt_y, color='blue', marker='.', linestyle='dashed')
        pyplot.title("Optimal Solution with cost: " + str(cost_opt))'''

        fig, axs = pyplot.subplots(1, 2, constrained_layout=True)

        axs[0].plot(x, y, color='red', marker='.', linestyle='dashed')
        axs[0].set_title("Solution with cost: " + str(cost))

        axs[1].plot(opt_x, opt_y, color='blue', marker='.', linestyle='dashed')
        axs[1].set_title("Optimal Solution with cost: " + str(cost_opt))

        pyplot.show()
        




def args_check() -> bool:
    if len(sys.argv) != 2+1 and len(sys.argv) != 3+1:
        print("Usage: python TspPlot <SOURCE_TSP_FILE> <SOLUTION_TOUR_FILE> (OPTIONAL)<OPTIMAL_SOLUTION_TOUR_FILE>")
        exit()

    for i in range(1, len(sys.argv), 1):
            if not os.path.isfile(sys.argv[i]):
                print("File \"" + sys.argv[i] + "\" does not seems to exist")
                exit()
    
    if len(sys.argv) == 2+1:
        return False
    else:
        return True
        


def get_points_from_tsp_file(filename:str) -> tuple[np.ndarray, np.ndarray, int]:
    with open(filename) as f:
        lines = f.readlines()
    
    supported_edge_weight_type : list[str] = ["EUC_2D","MAN_2D","MAX_2D", "ATT"]
    edge_wgt_type = -1

    while not "NODE_COORD_SECTION" in lines[0]:
        if "DIMENSION" in lines[0]:
            dim = int(lines[0].split(" : ")[1])
        elif "EDGE_WEIGHT_TYPE" in lines[0]:
            for i in range(len(supported_edge_weight_type)):
                if supported_edge_weight_type[i] in lines[0]:
                    edge_wgt_type = i
        del lines[0]
    del lines[0] # move to first line with coordinates

    x = np.zeros(shape=[dim+1], dtype=float)
    y = np.zeros(shape=[dim+1], dtype=float)

    i = 0
    while not "EOF" in lines[i]:
        nums = lines[i].split(' ')
        index, x[i], y[i] = [int(nums[0]), float(nums[1]), float(nums[2])]
        i += 1

    # close the tour
    x[x.shape[0]-1] = x[0]
    y[y.shape[0]-1] = y[0]

    return (x,y, edge_wgt_type)

def get_solution_from_tour_file(filename:str) -> np.ndarray:
    with open(filename) as f:
        lines = f.readlines()
    
    while not "TOUR_SECTION" in lines[0]:
        if "DIMENSION" in lines[0]:
            dim = int(lines[0].split(" : ")[1])
        del lines[0]
    del lines[0] # move to first line with coordinates
    
    seq = np.zeros(shape=[dim+1], dtype=int)

    i = 0
    while not "EOF" in lines[i] and not "-" in lines[i]:
        seq[i] = int(lines[i])
        i += 1
    
    seq -= 1

    # close the tour
    seq[seq.shape[0]-1] = seq[0]

    return seq

# arrays must be sorted according to solution
def compute_total_cost(x:np.ndarray, y:np.ndarray, edge_wgt_type:int) -> float:
    cost:float = 0.0

    if edge_wgt_type == 0:
        for i in range(x.shape[0]-1):
            cost += euc_2D(x[i],y[i], x[i+1],y[i+1])
    elif edge_wgt_type == 1:
        for i in range(x.shape[0]-1):
            cost += man_2D(x[i],y[i], x[i+1],y[i+1])
    elif edge_wgt_type == 2:
        for i in range(x.shape[0]-1):
            cost += max_2D(x[i],y[i], x[i+1],y[i+1])
    elif edge_wgt_type == 3:
        for i in range(x.shape[0]-1):
            cost += att(x[i],y[i], x[i+1],y[i+1])
    return cost

def euc_2D(x1:float, y1:float, x2:float, y2:float) -> float:
    return math.floor(math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)))

def man_2D(x1:float, y1:float, x2:float, y2:float) -> float:
    return math.floor(abs(x1-x2) + abs(y1-y2))

def max_2D(x1:float, y1:float, x2:float, y2:float) -> float:
    return math.floor(max(abs(x1-x2), abs(y1-y2)))

def att(x1:float, y1:float, x2:float, y2:float) -> float:
    return math.ceil(math.sqrt( ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)) / 10.0 ))

if __name__ == "__main__" :
    main()