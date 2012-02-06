# coding: utf-8
#
#  Copyright (C) 2008-2010 Istituto per l'Interscambio Scientifico I.S.I.
#  You can contact us by email (isi@isi.it) or write to:
#  ISI Foundation, Viale S. Severo 65, 10133 Torino, Italy. 
#
#  This program was written by André Panisson <panisson@gmail.com>
#
from pymobility.models.mobility import gauss_markov
'''
Created on Jan 24, 2012

@author: André Panisson
@contact: panisson@gmail.com
@organization: ISI Foundation, Torino, Italy
'''
import numpy as np
from models import truncated_levy_walk, random_direction, random_waypoint, random_walk
import logging

logging.basicConfig(format='%(asctime)-15s - %(message)s', level=logging.INFO)
logger = logging.getLogger("simulation")

DRAW = True
CALCULATE_CONTACTS = True

nr_nodes = 100
MAX_X, MAX_Y = 100, 100

#MIN_V, MAX_V = 0.003, 0.03
MIN_V, MAX_V = 0.1, 1.

MAX_WT = 100.
RANGE = 1.

STEPS_TO_IGNORE = 0#10000

trace_file = file("/tmp/trace.txt", 'w')

if DRAW:
    #import matplotlib
    #matplotlib.use('PS')
    import matplotlib.pyplot as plt
    plt.ion()
    ax = plt.subplot(111)
    line, = ax.plot(range(MAX_X), range(MAX_X), linestyle='', marker='.')

contacts = [None]*nr_nodes

def contacts_list():
    l = []
    for i in range(nr_nodes):
        if len(contacts[i])>0:
            for n in contacts[i]:
                l.append((i,n))
    return l

if DRAW and CALCULATE_CONTACTS:
    for l in range(100):
        ax.plot([], [], 'b-')
        
step = 0

# UNCOMMENT THE MODEL YOU WANT TO USE
#for x,y in truncated_levy_walk(nr_nodes, dimensions=(MAX_X, MAX_Y)):
#for x,y in random_direction(nr_nodes, dimensions=(MAX_X, MAX_Y)):
#for x,y in random_waypoint(nr_nodes, dimensions=(MAX_X, MAX_Y), velocity=(MIN_V, MAX_V), wt_max=MAX_WT):
#for x,y in random_walk(nr_nodes, dimensions=(MAX_X, MAX_Y)):
for x,y in gauss_markov(nr_nodes, dimensions=(MAX_X, MAX_Y), alpha=0.99):
    
    step += 1
    if step%10000==0: logger.info('Step %s'% step)
    if step < STEPS_TO_IGNORE: continue
    
    if CALCULATE_CONTACTS:
        
        for i in xrange(nr_nodes):
            d = np.sqrt(np.square(y-y[i]) + np.square(x-x[i]))
            contacts[i] = [n for n in np.where(d<RANGE)[0] if n!=i]
            
        trace_file.write(str(contacts_list())+"\n")
    
    if DRAW:
        
        if CALCULATE_CONTACTS:
            lnr = 1
            for i in xrange(nr_nodes):
                for j in contacts[i]:
                    if j > i:
                        ax.lines[lnr].set_data([x[i],x[j]], [y[i],y[j]])
                        lnr += 1
            for i in xrange(lnr, 100):
                ax.lines[i].set_data([],[])
        
        line.set_data(x, y)
        plt.draw()
    