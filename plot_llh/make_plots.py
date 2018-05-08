#!/usr/bin/python

import sys
sys.path.append('/Users/Teps/Theory')
#import header as h
#sys.path.append('/Users/Teps/Theory/HESE')
#import anarchy_header as ah
sys.path.append('/Users/Teps/Theory/HESE/Carlos')
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import MinimalTools as MT
import PhysConst as PC
import LVCPT as LVCPT
import numpy as np

import sys,os

rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})
cols = ['#29A2C6','#FF6D31','#FFCB18','#73B66B','#EF597B', '#333333']

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 18}

if len(sys.argv)< 2:
    print sys.argv
    print "Usage: make_plots.py input_filepath."
    exit(1)

#colors_schemes = ["Greens","Reds","Blues","PuRd","summer"]
colors_schemes = ["Greens","Reds","Blues","spring","summer"]
#colors_schemes = ["Greens","Reds","Blues","cool","summer"]
#colors_schemes = ["Greens","Reds","Blues","PuRd","summer"]
#colors_schemes = ["Blues","Greens","Reds","PuRd","summer"]
#colors_schemes = ["Greys","Greens","Reds","PuRd","summer"]
#colors_schemes = ["Greys","Greys","Greys","Greys","summer"]
#colors_schemes = ["PuRd","summer"]
output_format = "pdf"

# if True then will plot all files in the same figure
use_same_canvas = True 
figure = None

for i in range(len(sys.argv)-1):
    infile = sys.argv[i+1]
    print "Load data: " + str(infile)
    if infile[-3:] == 'txt':
        flavor_list = np.genfromtxt(infile)
    else:
        flavor_list = np.load(infile)
    if len(flavor_list[~np.isfinite(flavor_list)]) != 0:
        fl = []
        for x in flavor_list:
            if np.sum(~np.isfinite(x)) == 0:
                fl.append(x.tolist())
        flavor_list = np.array(fl)
    print flavor_list
    print "Done loading data"

    if not use_same_canvas :
        filename = "triangle_plot_"+os.path.basename(sys.argv[i+1])[:-4]
    else :
        filename = "triangle_plot"+os.path.basename(sys.argv[i+1])[:-4]

    # plots scale and diviions
    scale = 8
    divisions = 40

    print "Begin making plot ..."
    triangle_collection = LVCPT.GenerateTriangles(scale,divisions*2)
    map(lambda f : LVCPT.AddPointToTriangleCollection(f,triangle_collection, scale = scale, divisions = divisions*2),flavor_list)

    if use_same_canvas:
        figure = LVCPT.MakeFlavorTriangle(flavor_list, divisions = 5, save_file=True, 
                       filename = filename + "_hist", icolor = "g", icolormap = colors_schemes[i],
                       triangle_collection=triangle_collection, figure = figure, alpha = 0.8,
                       initial_flavor_ratio = [0,1,0], subdivisions = True, color_scale = "log",
                       output_format = output_format, inner_line_color ="silver",add_default_text = False,
                       plot_color_bar =True)

    else:
        figure = LVCPT.MakeFlavorTriangle(flavor_list, divisions = 5, save_file=True,
                       filename = filename + "_hist", icolor = "g", icolormap = colors_schemes[i],
                       triangle_collection=triangle_collection, alpha = 0.8,
                       initial_flavor_ratio = [0,1,0], subdivisions = True, color_scale = "log",
                       output_format = output_format, inner_line_color = "silver",add_default_text = False,
                       plot_color_bar =True)

    print "Done making plot: " + filename + "_hist."+output_format

ax = figure.get_axes()[0]
#ax.plot(6.5-0.35,5.5+0.3+0.35,"o", color = "grey", ms  = 10, markeredgewidth = 0.1, alpha = 0.9)
#ax.text(6.7-0.35,5.44+0.3+0.35,r"$(1-x:x:0)$", fontsize = 16)
#ax.axvline(x = 7.9)
fsz = 32
ax.plot(6.5-0.1,5.5+0.3+0.35+0.2+0.2,"o", color = "gold", ms  = 25, markeredgewidth = 0.1, alpha = 0.9)
ax.text(6.7-0.1,5.44+0.3+0.35+0.2+0.2,r"$(1:2:0)$", fontsize = fsz)

ax.plot(6.5-0.1,5.5+0.35+0.2,"o", color = "#2B653E", ms  = 25, markeredgewidth = 0.1, alpha = 0.9)
ax.text(6.7-0.1,5.44+0.35+0.2,r"$(1:0:0)$", fontsize = fsz)

ax.plot(6.5-0.1,5.5-0.3+0.35-0.2+0.2,"o", color = "#CA323D", ms  = 25, markeredgewidth = 0.1, alpha = 0.9)
ax.text(6.7-0.1,5.44-0.3+0.35-0.2+0.2,r"$(0:1:0)$", fontsize = fsz)

ax.plot(6.5-0.1,5.5-0.3+0.35-0.3-0.4+0.2,"o", color = "#2D4676", ms  = 25, markeredgewidth = 0.1, alpha = 0.9)
ax.text(6.7-0.1,5.44-0.3+0.35-0.3-0.4+0.2,r"$(0:0:1)$", fontsize = fsz)

plt.savefig("./plots/"+filename+"."+output_format, dpi = 600, bbox_inches='tight')

exit(1)

##os.system("cd plots")
##os.system("gs triangle_plot_hist.eps")
##os.system("cd ..")

