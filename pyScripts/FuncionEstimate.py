import numpy as np
import os
#import matplotlib
from matplotlib import pyplot as plt
import re
from scipy import signal
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

files_dir = "../run/nn_hyperparam_tuning_data"
cost_csv_fname = os.path.join(files_dir, "nn_almostbest_costs.csv")
iters_csv_fname = os.path.join(files_dir, "nn_almostbest_iterCount.csv")
nCores:int = 10


def read_csv(fname:str) -> [np.ndarray, np.array, np.array, np.ndarray]:
    with open(fname, 'r') as csvfile:
        firstline = csvfile.readline().strip().split(' ')
        ncols = int(firstline[0])
        hyperparam_vals = []
        for i in range(ncols):
            hyperparam_vals.append(float(re.findall(r'\d.\d+', firstline[i+1])[0]))
        nodes_num_list = []
        rows = []
        runtimes = []
        for row in csvfile:
            row = row.strip().split(' ')
            nodes_num_list.append(float(re.findall(r'\d+', row[0])[0]))
            rdata = []
            for j in range(ncols):
                rdata.append(float(row[j + 1]))
            rows.append(rdata)
            runtimes.append(float(row[ncols+1]))
    runtimes = np.array(runtimes)
    nodes_num_list =np.array(nodes_num_list)
    hyperparam_vals = np.array(hyperparam_vals)
    rows = np.array(rows)
    return (rows, nodes_num_list, hyperparam_vals, runtimes)

cost_data, nNodes, param_values, runtime_data = read_csv(cost_csv_fname)
iter_data, nNodes, param_values, runtime_data = read_csv(iters_csv_fname)

tot_iters = np.empty(shape=iter_data.shape, dtype=float)
for i in range(param_values.shape[0]):
    tot_iters[:,i] = iter_data[:,i] * runtime_data
np.round(tot_iters, out=tot_iters)

best_cost_pos = np.argmin(cost_data, axis=1)
best_totiters_data = tot_iters[np.arange(cost_data.shape[0]),best_cost_pos]
best_param_values = param_values[best_cost_pos]



ratio = best_totiters_data / nNodes
kernel_ratio = np.power(ratio, 1/3)



poly = PolynomialFeatures(degree=1, include_bias=False)
poly_features = poly.fit_transform(kernel_ratio.reshape(-1,1))
model = LinearRegression()
model.fit(poly_features, best_param_values)
print('Model coeffs\n' + str(model.coef_))
predicted = model.predict(poly_features)
print('predictions:\n' + str(predicted))
print('Avg Pred/True ratio = '+ str(np.average(np.abs(predicted / best_param_values))))

fig = plt.figure()

ax = plt.gca()

plt.scatter(ratio, predicted, c='blue', marker='.', s=16)
plt.scatter(ratio, best_param_values, c='red', marker='.', s=16)

plt.xlabel('number of iterations / number of nodes')
plt.ylabel('Grasp Chance')
# ax.grid(which='major')
# ax.grid(which='minor', linestyle='--')
plt.xscale('log')

plt.show()