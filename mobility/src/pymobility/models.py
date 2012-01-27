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
		
#		if i%100 == 0:
#			print np.average(velocity)

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

def random_contact(nr_nodes):
	a = np.array(range(nr_nodes))
	while True:
		i = np.random.randint(nr_nodes)
		j = i
		while j==i:
			j = np.random.randint(nr_nodes)
		yield a[i], a[j]
