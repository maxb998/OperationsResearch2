#!/usr/bin/env python2

from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import re
from scipy import signal

from optparse import OptionParser

# parameters
defLW = 1.2  # default line width
defMS = 2 # default marker size
dashes = ['-', '--', '-.', ':',	'-', '--', '-', '--', '-.', ':', '-', '--', '-', '--', '-.', ':', '-', '--', ':', '-.']
markers = ['+', 'x', '.', '^', 'o', 'd', 'p', '1', '.', 'x', '^', '+', 'x', '.', '^', 'o', 'd', 'p', '1', '.']
colors = ['r', 'b', 'y', 'g', 'm', 'c', 'orange', 'limegreen', 'midnightblue', 'violet', 'khaki', 'peru', 'pink', 'goldenrod', 'darkcyan', 'darkgreen', 'orangered', 'gray', 'darkorange', 'crimson']


class CmdLineParser(object):
	def __init__(self):
		self.parser = OptionParser(usage='usage: python2 perfprof.py [options] cvsfile.csv outputfile.pdf')
		# default options
		self.parser.add_option("-D", "--delimiter", dest="delimiter", default=';', help="delimiter for input files")
		self.parser.add_option("-M", "--maxratio", dest="maxratio", default=-1, type=float, help="maxratio for perf. profile")
		self.parser.add_option("-S", "--shift", dest="shift", default=0, type=float, help="shift for data")
		self.parser.add_option("-L", "--logplot", dest="logplot", action="store_true", default=False, help="log scale for x")
		self.parser.add_option("-T", "--timelimit", dest="timelimit", default=1e99, type=float, help="time limit for runs")
		self.parser.add_option("-P", "--plot-title", dest="plottitle", default=None, help="plot title")
		self.parser.add_option("-X", "--x-label", dest="xlabel", default='Ratio', help="x axis label")
		self.parser.add_option("-Y", "--y-label", dest="ylabel", default='Number of Nodes', help="y axis label")
		self.parser.add_option("-B", "--bw", dest="bw", action="store_true", default=False, help="plot B/W")
		self.parser.add_option("--sortNNodes", dest="sortNNodes", action="store_true", default=False, help="Sortes the graph y axis using the number of nodes of each instance")
		self.parser.add_option("--smooth", dest="smooth", action="store_true", default=False, help="Smooth the output")

	def addOption(self, *args, **kwargs):
		self.parser.add_option(*args, **kwargs)

	def parseArgs(self):
		(options, args) = self.parser.parse_args()
		options.input = args[0]
		options.output = args[1]
		return options


def readTable(fp, delimiter):
	"""
	read a CSV file with performance profile specification
	the format is as follows:
	ncols algo1 algo2 ...
	nome_istanza tempo(algo1) tempo(algo2) ...
	...
	"""
	firstline = fp.readline().strip().split(delimiter)
	ncols = int(firstline[0])
	assert(ncols <= len(markers))
	cnames = firstline[1:]
	rnames = []
	rows = []
	for row in fp:
		row = row.strip().split(delimiter)
		rnames.append(row[0])
		rdata = np.empty(ncols)
		for j in range(ncols):
			if row[j+1] == 'null':
				rdata[j] = np.Inf
			else:
				rdata[j] = float(row[j + 1])
		rows.append(rdata)
	data = np.array(rows)
	return (rnames, cnames, data)


def main():
	parser = CmdLineParser()
	opt = parser.parseArgs()
	print(opt)
	# read data
	rnames, cnames, data = readTable(open(opt.input, 'r'), opt.delimiter)
	nrows, ncols = data.shape
	# add shift
	data = data + opt.shift
	# compute ratios
	minima = data.min(axis=1)
	ratio = data
	for j in range(ncols):
		ratio[:, j] = data[:, j] / minima
	# compute maxratio
	if opt.maxratio == -1:
		opt.maxratio = ratio.max() * 1.01
	# any time >= timelimit will count as maxratio + bigM (so that it does not show up in plots)
	for i in range(nrows):
		for j in range(ncols):
			if data[i,j] >= opt.timelimit:
				ratio[i,j] = opt.maxratio + 1e6
    
	y = np.arange(nrows, dtype=np.float64) / nrows

	if opt.sortNNodes:
		for i in range(len(rnames)):
			y[i] = re.findall(r'\d+', rnames[i])[0]

		if opt.smooth:
			windowSizes = [25, 19, 13, 7, 3]
			polyOrders = [3, 3, 3, 2, 2]
			for i in range(ncols):
				cpRatio = np.copy(ratio[:,i])
				notInfiniteMap = cpRatio < np.Infinity
				for j in range(len(windowSizes)):
					if cpRatio[notInfiniteMap].shape[0] > windowSizes[j]:
						cpRatio[notInfiniteMap] = signal.savgol_filter(cpRatio[notInfiniteMap], windowSizes[j], polyOrders[j])
						break
				
				ratio[:,i] = cpRatio


	else:
		# sort ratios
		ratio.sort(axis=0)
	
	# plot first
	for j in range(ncols):
		options = dict(label=cnames[j],
				linewidth=defLW, linestyle=dashes[j],
				marker=markers[j], markeredgewidth=defLW, markersize=defMS)
		#plt.step(ratio[:,j], y, label=cnames[j], linewidth=defLW, marker=markers[j], markersize=defMS)
		if opt.bw:
			options['markerfacecolor'] = 'w'
			options['markeredgecolor'] = 'k'
			options['color'] = 'k'
		else:
			options['color'] = colors[j]
		if opt.sortNNodes:
			plt.semilogy(ratio[:, j], y, **options)
		elif opt.logplot:
			plt.semilogx(ratio[:, j], y, **options)
		else:
			plt.plot(ratio[:, j], y, **options)
	if not opt.sortNNodes:
		plt.axis([1, opt.maxratio, 0, 1])
	plt.legend(loc='lower right')
	if opt.plottitle is not None:
		plt.title(opt.plottitle)
	plt.xlabel(opt.xlabel)
	plt.ylabel(opt.ylabel)
	plt.savefig(opt.output)

if __name__ == '__main__':
	main()