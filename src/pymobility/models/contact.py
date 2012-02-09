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

import unittest
from numbers import Number

'''
Created on Jan 31, 2012

@author: André Panisson
@contact: panisson@gmail.com
@organization: ISI Foundation, Torino, Italy
'''
import numpy as np

def random_contact(nr_nodes):
    a = np.array(range(nr_nodes))
    while True:
        i = np.random.randint(nr_nodes)
        j = i
        while j==i:
            j = np.random.randint(nr_nodes)
        yield a[i], a[j]

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
    for i in range(len(partner)):
        if len(partner[i])>0:
            for k in partner[i]:
                l.append((i,k))
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
        
class Test(unittest.TestCase):
    
    def test_modelB(self):
        
        N=100                  # number of agents (going from 0 to N-1)              
        b0 = .51               # transition 0 -> 1 parameter       
        b1 = .86               # state change parameter if the agent is active
        mu = .95               # aggregating tendance parameter
        
        PI = lambda x,y : 1./(1+1.*(x-y)/N)
        f = lambda x,y : 1./(1+1.*(x-y)/N)
        
        t = 0
        for contacts in modelB(b0,b1,mu,PI,f,N,seed_value=0xbbbbbb):
            t += 1
            if t==10000:
                assert contacts == [(0, 91), (3, 16), (5, 98), (8, 84), (16, 3), (21, 54), (24, 47), (27, 78), (36, 72), (47, 24), (52, 94), (54, 21), (55, 61), (58, 59), (59, 58), (61, 55), (72, 36), (78, 27), (84, 8), (91, 0), (94, 52), (98, 5)]
                break
            
    def test_model_het(self):
        
        N=100                  # number of agents (going from 0 to N-1)
        mu = .86               # aggregating tendance parameter
        
        PI = lambda x,y : 1./(1+1.*(x-y)/N)
        f = lambda x,y : 1./(1+1.*(x-y)/N)
        f2 = lambda x : 1./(1+1.*(x)/N)
        
        t = 0
        for contacts in model_het(f2,mu,PI,f,N,seed_value=0xbbbbbb):
            t += 1
            if t==10000:
                assert contacts == [(0, 36), (1, 21), (2, 4), (3, 92), (4, 2), (5, 70), (5, 87), (6, 46), (7, 53), (8, 80), (8, 52), (9, 15), (10, 27), (11, 23), (12, 90), (13, 81), (14, 34), (15, 9), (18, 35), (19, 75), (20, 98), (21, 1), (22, 78), (23, 11), (24, 58), (25, 26), (26, 25), (27, 10), (28, 31), (29, 95), (30, 50), (31, 28), (32, 79), (33, 68), (34, 14), (35, 18), (36, 0), (37, 44), (38, 47), (39, 55), (40, 74), (41, 63), (42, 73), (43, 93), (44, 37), (45, 67), (46, 6), (47, 38), (48, 97), (49, 96), (50, 30), (51, 66), (52, 8), (52, 80), (53, 7), (54, 61), (55, 39), (57, 82), (58, 24), (60, 83), (61, 54), (62, 91), (63, 41), (65, 89), (66, 51), (67, 45), (68, 33), (70, 5), (70, 87), (73, 42), (74, 40), (75, 19), (76, 99), (77, 94), (78, 22), (79, 32), (80, 8), (80, 52), (81, 13), (82, 57), (83, 60), (87, 5), (87, 70), (89, 65), (90, 12), (91, 62), (92, 3), (93, 43), (94, 77), (95, 29), (96, 49), (97, 48), (98, 20), (99, 76)]
                break
