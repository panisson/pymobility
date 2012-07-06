# coding: utf-8
#
#  Copyright (C) 2012 Istituto per l'Interscambio Scientifico I.S.I.
#  You can contact us by email (isi@isi.it) or write to:
#  ISI Foundation, Viale S. Severo 65, 10133 Torino, Italy. 
#
#  This program was written by André Panisson <panisson@gmail.com>
#
####### modelB implementation #######
#
# Author : Juliette Stehle
# Reference publication: 
#       Stehle, J., Barrat, A. & Bianconi, G. 
#       Dynamical and bursty interactions in social networks. Physical Review E 81, 1-4 (2010).
# Last change : May 17 2010
#
######################################

from numbers import Number

'''
Created on Jan 31, 2012

@author: André Panisson
@contact: panisson@gmail.com
@organization: ISI Foundation, Torino, Italy
'''
import numpy as np

def dynamic_gnp(n, p):
    '''
    Implementation of the Dynamic G(n,p) graph as discussed in the following paper:
        Andrea E. F. Clementi, Francesco Pasquale, Angelo Monti, and Riccardo Silvestri. 2007. 
        Communication in dynamic radio networks. In Proceedings of the twenty-sixth annual 
        ACM symposium on Principles of distributed computing (PODC '07). ACM, New York, NY, USA, 205-214.
    At each time slot t of the execution of the protocol, a (new) graph G(t) is selected
    according to the well-known random graph model G(n,p) where n is the number of nodes and p is the
    edge probability.
    
    Required arguments:
    
      *n*:
        The number of vertices in the graph.

      *p*   
        The probability for drawing an edge between two arbitrary vertices (G(n,p) graph)
    
    '''
    while True:
        m = np.random.rand(n, n)
        c = np.where(m < p)
        yield [(i,j) for i,j in zip(c[0], c[1]) if i < j]

def dynamic_gnm(n, m):
    '''
    Implementation of the Dynamic G(n,m) graph as discussed in the following paper:
    
        Andrea E. F. Clementi, Francesco Pasquale, Angelo Monti, and Riccardo Silvestri. 2007. 
        Communication in dynamic radio networks. In Proceedings of the twenty-sixth annual 
        ACM symposium on Principles of distributed computing (PODC '07). ACM, New York, NY, USA, 205-214.
    
    At each time slot t of the execution of the protocol, a (new) graph G(t) is selected
    according to the well-known random graph model G(n,m) where n is the number of nodes and m is the
    number of edges.
    
    Required arguments:
    
      *n*:
        The number of vertices in the graph.

      *m*:    
        The number of edges in the graph (for G(n,m) graphs).
    
    '''
    while True:
        contacts = []
        while len(contacts) < m:
            i = np.random.randint(n)
            j = np.random.randint(n)
            if i == j or (i, j) in contacts: continue
            contacts.append((i, j))
        yield contacts
        
def edge_markovian(n, p, q, g=0):
    '''
    Implementation of the edge-Markovian dynamic graph as discussed in the following paper:
    
    Andrea E.F. Clementi, Claudio Macci, Angelo Monti, Francesco Pasquale, and Riccardo Silvestri. 2008. 
        Flooding time in edge-Markovian dynamic graphs. In Proceedings of the 
        twenty-seventh ACM symposium on Principles of distributed computing (PODC '08). 
        ACM, New York, NY, USA, 213-222.
    
    Starting from an arbitrary initial edge probability distribu- 
    tion, at every time step, every edge changes its state (exist- 
    ing or not) according to a two-state Markovian process with 
    probabilities p (edge birth-rate) and q (edge death-rate). If 
    an edge exists at time t then, at time t + 1, it dies with prob- 
    ability q. If instead the edge does not exist at time t, then 
    it will come into existence at time t + 1 with probability p. 
    
    Required arguments:
    
      *n*:
        The number of vertices in the graph.

      *p*
        If the edge does not exist at time t then it will come into existence at time
        t + 1 with probability p.
        
      *p*
        If an edge exists at time t then, at time t + 1, it dies with probability q.
    
      *g*:
        an arbitrary initial probability distribution over the set [n] yielding E0.

    
    '''
    # adjacency matrix
    a = np.zeros((n,n))
    
    # initial set of edges
    a[np.where(np.random.rand(n, n) < g)] = 1
    
    while True:
        m = np.random.rand(n, n)
        up = np.where(np.logical_and(a == 0, m < p))
        down = np.where(np.logical_and(a > 0, m < q))
        a[up] = 1.
        a[down] = 0.
        
        c = np.nonzero(a)
        yield [(i,j) for i,j in zip(c[0], c[1]) if i < j]
        
def continuous_time_edge_markovian(n, lmbd):
    '''
    Implementation of the continuous-time edge-Markovian dynamic graph as discussed in the following paper:
    
    Augustin Chaintreau, Abderrahmen Mtibaa, Laurent Massoulie, and Christophe Diot. 2007. 
        The diameter of opportunistic mobile networks. In Proceedings of the 
        2007 ACM CoNEXT conference (CoNEXT '07). ACM, New York, NY, USA, , Article 12 , 12 pages.
    
    We assume that, for any pairs of nodes (u, v), the times of
    contact are separated by exponential random variables.
    '''
    from pymobility.models.mobility import E
    a = np.zeros((n,n))
    a = E(lmbd, a)
    
    while True:
        a -= 1.
        c = np.where(a <= 0.)
        yield [(i,j) for i,j in zip(c[0], c[1]) if i < j]
        a[c] = E(lmbd, c[0])
        
def broad_continuous_time_edge_markovian(n, alpha):
    '''
    This model is similar to the continuous-time edge-Markovian dynamic graph as discussed in the following paper:
    
    Augustin Chaintreau, Abderrahmen Mtibaa, Laurent Massoulie, and Christophe Diot. 2007. 
        The diameter of opportunistic mobile networks. In Proceedings of the 
        2007 ACM CoNEXT conference (CoNEXT '07). ACM, New York, NY, USA, , Article 12 , 12 pages.
        
    The difference is in the generated inter-contact times.    
    We assume that, for any pairs of nodes (u, v), the times of
    contact are separated by random variables from a power-law distribution.
    '''
    # define a Power Law Distribution
    P = lambda ALPHA, SCALE, SAMPLES: (SAMPLES ** (1./(1.-alpha))) * SCALE
    
    a = P(alpha, 1., np.random.rand(n,n))
    
    while True:
        a -= 1.
        c = np.where(a <= 0.)
        yield [(i,j) for i,j in zip(c[0], c[1]) if i < j]
        a[c] = P(alpha, 1., np.random.rand(*c[0].shape))

def mobility_contact(mobility_model, contact_range=1.0):
    '''
    Contact model based on a mobility model.
    At each time step, this model gets the position of the nodes accordingly to the given mobility model
    and calculates the distance of each node pair.
    The list of contacts generated in each time step is the list of node pairs with distance lower than *contact_range*.
    
    Required arguments:
    
      *nr_nodes*:
        Integer, the number of nodes.
      
      *contact_range*:
        Double, the maximum euclidean distance in which two nodes are in contact.
    '''
    from scipy.spatial.distance import cdist
    for xy in mobility_model:
        d = cdist(xy,xy)
        contacts = []
        for (i,j) in zip(*np.where(d<contact_range)):
            if j > i: contacts.append((i,j))
        yield contacts

##################################### 
def __inactive_partner_choice(t_lc,t,inactive_agents,PI,eta): 
    from random import random
    """Returns an inactive partner with a probability 
    proportional to 1/(1+t-t[i]) where t[i] is the last state change time"""
    L=len(inactive_agents)
    proba,cumulated_prob,Z = [0]*L, [0]*(L+1),0               
    for j in range(0,L): 
        c = eta[j] if not isinstance(eta,Number) else 1.
        proba[j] = c*PI(t,t_lc[inactive_agents[j]]) 
        Z = Z + proba[j] 
    for j in range(1,L+1): 
        cumulated_prob[j]= cumulated_prob[j-1] + proba[j-1]/Z 
    k = random() 
    j = 0 
    while k > cumulated_prob[j+1]: j += 1 
    return inactive_agents[j] 

#####################################

##################################### 
def __update_state(i,b1,mu,PI,f,t,t_lc,partner,inactive_agents): 
    """Run the change of an active agent : if the state changes,
    it can whether leave the group, whether introduce an inactive agent."""
    from random import random

    k1 = random() 
    p = len(partner[i])
    ti = t_lc[i]
    
    
    if isinstance(b1,Number):
        eta = b1 # modelb
    else:
        eta = 1-b1[i] # model_het
    
    if k1 <= eta*f(t,ti): 
        k2 = random()
 
        ### The agent i leaves the group 
        if k2 <= mu:

            for key in partner[i]:
                t_lc[key] = t 
                del partner[key][i] 
            if p==1 : inactive_agents.extend(partner[i])
            inactive_agents.append(i)
            t_lc[i],partner[i] = t,{} 

        ### The agent i introduces an inactive agent j to the group 
        elif len(inactive_agents)> 1:
            j = __inactive_partner_choice(t_lc,t,inactive_agents,PI,b1)
            inactive_agents.remove(j)
            partner[j] = partner[i].copy() 
            partner[j][i]=t 
            for key in partner[j]: 
                partner[key][j]=t 
            t_lc[j] = t 
            for key in partner[j]: 
                t_lc[key] = t
#####################################

def __contacts_list(partner):
    l = []
    for i in range(len(partner)): l.extend((i,k) for k in partner[i])
    return l

##################################### 
def modelB(b0, b1, mu, PI, f, N, seed_value=None):
    from random import seed,randint,random
    """Realize an execution of the model for a parameter set"""

    ### Variables initialization 
    last_change_time = [0]*N    # Time of the last state change of each agent 
    partner = [{} for _ in range(N)]  # List of interacting partners of each agent 
                                #         with the time of the beginning of the contact 
    inactive_agents = range(0,N)

    seed(seed_value)
    t = 0 
    while True:
        i = randint(0,N-1) 
        if len(partner[i]) == 0 : # i is not interacting
            k = random()
            ti = last_change_time[i]
            if k <= b0*f(t,ti) and len(inactive_agents)!=1 :
                inactive_agents.remove(i)
                j = __inactive_partner_choice(last_change_time,t,inactive_agents,PI,b1) 
                partner[i]={j:t} 
                partner[j]={i:t}
                inactive_agents.remove(j)
                last_change_time[i] = t 
                last_change_time[j] = t
                
        else:                     # i is interacting
            __update_state(i,b1,mu,PI,f,t,last_change_time,
                         partner,inactive_agents)
        t += 1
        
        yield __contacts_list(partner)
        
##################################### 
def model_het(f2,mu,PI,f,N,seed_value=None):# PARAM,N,tmax):
    from random import seed,randint,random
    """Realize an execution of the model for a parameter set"""

    ### Variables initialization 
    last_change_time = [0]*N    # Time of the last state change of each agent 
    partner = [{} for _ in range(N)]   # List of interacting partners of each agent 
                                #         with the time of the beginning of the contact 
    inactive_agents = range(0,N)
    eta = [0]*N
    for i in range(N):
        eta[i] = f2(i)

    seed(seed_value)  
    t = 0 
    while True: 
        i = randint(0,N-1) 
        if len(partner[i]) == 0 : # i is not interacting
            k = random()
            ti = last_change_time[i]
            if k <= eta[i]*f(t,ti) and len(inactive_agents)!=1 :
                inactive_agents.remove(i)
                j = __inactive_partner_choice(last_change_time,t,inactive_agents,PI,eta) 
                partner[i]={j:t} 
                partner[j]={i:t}
                inactive_agents.remove(j)
                last_change_time[i] = t 
                last_change_time[j] = t
                
        else:                     # i is interacting
            __update_state(i,eta,mu,PI,f,t,last_change_time,
                         partner,inactive_agents)
        t += 1
        
        yield __contacts_list(partner)
        
