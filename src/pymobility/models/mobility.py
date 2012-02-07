# coding: utf-8
#
#  Copyright (C) 2008-2010 Istituto per l'Interscambio Scientifico I.S.I.
#  You can contact us by email (isi@isi.it) or write to:
#  ISI Foundation, Viale S. Severo 65, 10133 Torino, Italy. 
#
#  This program was written by André Panisson <panisson@gmail.com>
#
'''
Created on Jan 24, 2012

@author: André Panisson
@contact: panisson@gmail.com
@organization: ISI Foundation, Torino, Italy
'''
import numpy as np
from numpy.random import rand

# define a Uniform Distribution
U = lambda MIN, MAX, SAMPLES: rand(len(SAMPLES)) * (MAX - MIN) + MIN

# define a Truncated Power Law Distribution
P = lambda ALPHA, MIN, MAX, SAMPLES: ((MAX ** (ALPHA+1.) - 1.) * rand(len(SAMPLES)) + 1.) ** (1./(ALPHA+1.))

def random_waypoint(nr_nodes, dimensions, velocity=(0.1, 1.), wt_max=None):
    '''
    Random Waypoint model.
    
    Required arguments:
    
      *nr_nodes*:
        Integer, the number of nodes.
      
      *dimensions*:
        Tuple of Integers, the x and y dimensions of the simulation area.
      
    keyword arguments:
    
      *velocity*:
        Tuple of Integers, the minimum and maximum values for node velocity.
      
      *wt_max*:
        Integer, the maximum wait time for node pauses.
        If wt_max is 0 or None, there is no pause time.
    '''
    
    MAX_X,MAX_Y = dimensions
    MIN_V, MAX_V = velocity
    
    NODES = np.arange(nr_nodes)
    x = U(0, MAX_X, NODES)
    y = U(0, MAX_Y, NODES)
    x_waypoint = U(0, MAX_X, NODES)
    y_waypoint = U(0, MAX_Y, NODES)
    wt = np.zeros(nr_nodes)
    velocity = U(MIN_V, MAX_V, NODES)
        
    while True:
        # update node position
        theta = np.arctan2(y_waypoint - y, x_waypoint - x)
        x = x + velocity * np.cos(theta)
        y = y + velocity * np.sin(theta)
        # calculate distance to waypoint
        d = np.sqrt(np.square(y_waypoint-y) + np.square(x_waypoint-x))
        # update info for arrived nodes
        arrived = np.where(d<=velocity)[0]
        
        if wt_max:
            velocity[arrived] = 0.
            wt[arrived] = U(0, wt_max, arrived)
            # update info for paused nodes
            wt[np.where(velocity==0.)[0]] -= 1.
            # update info for moving nodes
            arrived = np.where(np.logical_and(velocity==0., wt<0.))[0]
        
        x_waypoint[arrived] = U(0, MAX_X, arrived)
        y_waypoint[arrived] = U(0, MAX_Y, arrived)
        velocity[arrived] = U(MIN_V, MAX_V, arrived)
        
        yield x,y
        
#        if i%100 == 0:
#            print np.average(velocity)

def random_walk(nr_nodes, dimensions, velocity=1., distance=1.):
    if velocity>distance:
        # In this implementation, each step is 1 second,
        # it is not possible to have a velocity larger than the distance
        raise Exception('Velocity must be <= Distance')
    
    fl = np.zeros(nr_nodes)+distance
    vel = np.zeros(nr_nodes)+velocity
    
    FL_DISTR = lambda SAMPLES: np.array(fl[:len(SAMPLES)])
    VELOCITY_DISTR = lambda FD: np.array(vel[:len(FD)])
    
    return stochastic_walk(nr_nodes, dimensions, FL_DISTR, VELOCITY_DISTR)

def truncated_levy_walk(nr_nodes, dimensions, FL_EXP=-2.6, FL_MAX=50., WT_EXP=-1.8, WT_MAX=100.):
    
    FL_DISTR = lambda SAMPLES: P(FL_EXP, 1., FL_MAX, SAMPLES)
    if WT_EXP and WT_MAX:
        WT_DISTR = lambda SAMPLES: P(WT_EXP, 1., WT_MAX, SAMPLES)
    else:
        WT_DISTR = None
    VELOCITY_DISTR = lambda FD: np.sqrt(FD)/10.
    
    return stochastic_walk(nr_nodes, dimensions, FL_DISTR, VELOCITY_DISTR, WT_DISTR)

def heterogeneous_truncated_levy_walk(nr_nodes, dimensions, WT_EXP=-1.8, WT_MAX=100., FL_EXP=-2.6, FL_MAX=50.):
    
    NODES = np.arange(nr_nodes)
    FL_MAX = P(-1.8, 10., FL_MAX, NODES)
    FL_MIN = FL_MAX/10.
    
    FL_DISTR = lambda SAMPLES: rand(len(SAMPLES)) * (FL_MAX[SAMPLES] - FL_MIN[SAMPLES]) + FL_MIN[SAMPLES]
    WT_DISTR = lambda SAMPLES: P(WT_EXP, 1., WT_MAX, SAMPLES)
    VELOCITY_DISTR = lambda FD: np.sqrt(FD)/10.
    
    return stochastic_walk(nr_nodes, dimensions, FL_DISTR, VELOCITY_DISTR, WT_DISTR)

def random_direction(nr_nodes, dimensions, wt_max=None, velocity=(0.1, 1.)):
    
    MIN_V, MAX_V = velocity
    FL_MAX = max(dimensions)
    
    FL_DISTR = lambda SAMPLES: U(0, FL_MAX, SAMPLES)
    if wt_max:
        WT_DISTR = lambda SAMPLES: U(0, wt_max, SAMPLES)
    else:
        WT_DISTR = None
    VELOCITY_DISTR = lambda FD: U(MIN_V, MAX_V, FD)
    
    return stochastic_walk(nr_nodes, dimensions, FL_DISTR, VELOCITY_DISTR, WT_DISTR)

def stochastic_walk(nr_nodes, dimensions, FL_DISTR, VELOCITY_DISTR, WT_DISTR=None):
    '''
    Base implementation for models with direction uniformly chosen from [0,pi]:
    random_direction, random_walk, truncated_levy_walk
    
    Required arguments:
    
      *nr_nodes*:
        Integer, the number of nodes.
      
      *dimensions*:
        Tuple of Integers, the x and y dimensions of the simulation area.
        
      *FL_DISTR*:
        A function that, given a set of samples, 
         returns another set with the same size of the input set.
        This function should implement the distribution of flight lenghts
         to be used in the model.
         
      *VELOCITY_DISTR*:
        A function that, given a set of flight lenghts, 
         returns another set with the same size of the input set.
        This function should implement the distribution of velocities
         to be used in the model, as random or as a function of the flight lenghts.
      
    keyword arguments:
    
      *WT_DISTR*:
        A function that, given a set of samples, 
         returns another set with the same size of the input set.
        This function should implement the distribution of wait times
         to be used in the node pause.
        If WT_DISTR is 0 or None, there is no pause time.
    '''
    
    MAX_X, MAX_Y = dimensions
    NODES = np.arange(nr_nodes)
    x = U(0, MAX_X, NODES)
    y = U(0, MAX_Y, NODES)
    fl = FL_DISTR(NODES)
    velocity = VELOCITY_DISTR(fl)
    theta = U(0, 2*np.pi, NODES)
    costheta = np.cos(theta)
    sintheta = np.sin(theta)
    wt = np.zeros(nr_nodes)
    
    while True:

        x = x + velocity * costheta
        y = y + velocity * sintheta
        
        # node bounces on the margins
        b = np.where(x<0)[0]
        x[b] = - x[b]; costheta[b] = -costheta[b]
        b = np.where(x>MAX_X)[0]
        x[b] = 2*MAX_X - x[b]; costheta[b] = -costheta[b]
        b = np.where(y<0)[0]
        y[b] = - y[b]; sintheta[b] = -sintheta[b]
        b = np.where(y>MAX_Y)[0]
        y[b] = 2*MAX_Y - y[b]; sintheta[b] = -sintheta[b]

        # update info for arrived nodes
        fl = fl - velocity
        arrived = np.where(np.logical_and(velocity>0., fl<=0.))[0]
        
        if WT_DISTR:
            velocity[arrived] = 0.
            wt[arrived] = WT_DISTR(arrived)
            # update info for paused nodes
            wt[np.where(velocity==0.)[0]] -= 1.
            # update info for moving nodes
            arrived = np.where(np.logical_and(velocity==0., wt<0.))[0]
            
        theta = U(0, 2*np.pi, arrived)
        costheta[arrived] = np.cos(theta)
        sintheta[arrived] = np.sin(theta)
        fl[arrived] = FL_DISTR(arrived)
        velocity[arrived] = VELOCITY_DISTR(fl[arrived])

        yield x,y


def gauss_markov(nr_nodes, dimensions, velocity_mean=1., alpha=1., variance=1.):
    '''
    Gauss-Markov Mobility Model, as proposed in 
    Camp, T., Boleng, J. & Davies, V. A survey of mobility models for ad hoc network research. 
    Wireless Communications and Mobile Computing 2, 483-502 (2002).
    
    Required arguments:
    
      *nr_nodes*:
        Integer, the number of nodes.
      
      *dimensions*:
        Tuple of Integers, the x and y dimensions of the simulation area.
        
    keyword arguments:
    
      *velocity_mean*:
        The mean velocity
        
      *alpha*:
        The tuning parameter used to vary the randomness
        
      *variance*:
        The randomness variance
    '''
    
    MAX_X, MAX_Y = dimensions
    NODES = np.arange(nr_nodes)
    x = U(0, MAX_X, NODES)
    y = U(0, MAX_Y, NODES)
    velocity =  np.zeros(nr_nodes)+velocity_mean
    theta = U(0, 2*np.pi, NODES)
    angle_mean = theta
    
    alpha2 = 1.0 - alpha
    alpha3 = np.sqrt(1.0 - alpha * alpha) * variance
    
    while True:

        x = x + velocity * np.cos(theta)
        y = y + velocity * np.sin(theta)
        
        # node bounces on the margins
        b = np.where(x<0)[0]
        x[b] = - x[b]; theta[b] = np.pi-theta[b]; angle_mean[b] = np.pi-angle_mean[b]
        b = np.where(x>MAX_X)[0]
        x[b] = 2*MAX_X - x[b]; theta[b] = np.pi-theta[b]; angle_mean[b] = np.pi-angle_mean[b]
        b = np.where(y<0)[0]
        y[b] = - y[b]; theta[b] = -theta[b]; angle_mean[b] = -angle_mean[b]
        b = np.where(y>MAX_Y)[0]
        y[b] = 2*MAX_Y - y[b]; theta[b] = -theta[b]; angle_mean[b] = -angle_mean[b]
        
        # calculate new speed and direction based on the model
        velocity = (alpha * velocity +
                    alpha2 * velocity_mean +
                    alpha3 * np.random.normal(0.0, 1.0, nr_nodes))
    
        theta = (alpha * theta +
                    alpha2 * angle_mean +
                    alpha3 * np.random.normal(0.0, 1.0, nr_nodes))
        
        yield x,y
        
def reference_point_group(nr_nodes, dimensions, velocity=(0.1, 1.), aggregation=0.1):
    '''
    Reference Point Group Mobility model.
    In this implementation, group trajectories follow a random direction model,
    while nodes follow a random walk around the group center.
    The parameter 'aggregation' controls how close the nodes are to the group center.
    
    Required arguments:
    
      *nr_nodes*:
        list of integers, the number of nodes in each group.
      
      *dimensions*:
        Tuple of Integers, the x and y dimensions of the simulation area.
        
    keyword arguments:
    
      *velocity*:
        Tuple of Doubles, the minimum and maximum values for group velocity.
        
      *aggregation*:
        Double, parameter (between 0 and 1) used to aggregate the nodes in the group.
        Usually between 0 and 1, the more this value approximates to 1,
        the nodes will be more aggregated and closer to the group center.
        With a value of 0, the nodes are randomly distributed in the simulation area.
        With a value of 1, the nodes are close to the group center.
    '''
    
    try:
        iter(nr_nodes)
    except TypeError:
        nr_nodes = [nr_nodes]
    
    NODES = np.arange(sum(nr_nodes))
    
    groups = []
    prev = 0
    for (i,n) in enumerate(nr_nodes):
        groups.append(np.arange(prev,n+prev))
        prev += n
    
    g_ref = np.empty(sum(nr_nodes), dtype=np.int)
    for (i,g) in enumerate(groups):
        for n in g:
            g_ref[n] = i
    
    FL_MAX = max(dimensions)
    MIN_V,MAX_V = velocity
    FL_DISTR = lambda SAMPLES: U(0, FL_MAX, SAMPLES)
    VELOCITY_DISTR = lambda FD: U(MIN_V, MAX_V, FD)
    
    MAX_X, MAX_Y = dimensions
    x = U(0, MAX_X, NODES)
    y = U(0, MAX_Y, NODES)
    fl = np.ones(sum(nr_nodes))
    velocity = VELOCITY_DISTR(fl)
    theta = U(0, 2*np.pi, NODES)
    costheta = np.cos(theta)
    sintheta = np.sin(theta)
    
    GROUPS = np.arange(len(groups))
    g_fl = FL_DISTR(GROUPS)
    g_velocity = VELOCITY_DISTR(g_fl)
    g_theta = U(0, 2*np.pi, GROUPS)
    g_costheta = np.cos(theta)
    g_sintheta = np.sin(theta)
        
    while True:

        x = x + velocity * costheta
        y = y + velocity * sintheta
        
        for (i,g) in enumerate(groups):
            
            # step to group direction + step to group center
            y_g = y[g]
            x_g = x[g]
            c_theta = np.arctan2(np.average(y_g) - y_g, np.average(x_g) - x_g)
            
            x[g] = x_g + g_velocity[i] * g_costheta[i] + aggregation*np.cos(c_theta)
            y[g] = y_g + g_velocity[i] * g_sintheta[i] + aggregation*np.sin(c_theta)
            
        # node and group bounces on the margins
        b = np.where(x<0)[0]
        x[b] = - x[b]; costheta[b] = -costheta[b]
        g_idx = np.unique(g_ref[b]); g_costheta[g_idx] = -g_costheta[g_idx]
        b = np.where(x>MAX_X)[0]
        x[b] = 2*MAX_X - x[b]; costheta[b] = -costheta[b]
        g_idx = np.unique(g_ref[b]); g_costheta[g_idx] = -g_costheta[g_idx]
        b = np.where(y<0)[0]
        y[b] = - y[b]; sintheta[b] = -sintheta[b]
        g_idx = np.unique(g_ref[b]); g_sintheta[g_idx] = -g_sintheta[g_idx]
        b = np.where(y>MAX_Y)[0]
        y[b] = 2*MAX_Y - y[b]; sintheta[b] = -sintheta[b]
        g_idx = np.unique(g_ref[b]); g_sintheta[g_idx] = -g_sintheta[g_idx]

        # update info for arrived nodes
        fl = fl - velocity
        arrived = np.where(np.logical_and(velocity>0., fl<=0.))[0]
        
        theta = U(0, 2*np.pi, arrived)
        costheta[arrived] = np.cos(theta)
        sintheta[arrived] = np.sin(theta)
        fl[arrived] = np.ones(len(arrived))
        velocity[arrived] = VELOCITY_DISTR(fl[arrived])
        
        # update info for arrived groups
        g_fl = g_fl - g_velocity
        g_arrived = np.where(np.logical_and(g_velocity>0., g_fl<=0.))[0]
        
        g_theta = U(0, 2*np.pi, g_arrived)
        g_costheta[g_arrived] = np.cos(g_theta)
        g_sintheta[g_arrived] = np.sin(g_theta)
        g_fl[g_arrived] = FL_DISTR(g_arrived)
        g_velocity[g_arrived] = VELOCITY_DISTR(fl[g_arrived])

        yield x,y
