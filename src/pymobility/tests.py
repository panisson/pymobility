# coding: utf-8
#
#  Copyright (C) 2012 Istituto per l'Interscambio Scientifico I.S.I.
#  You can contact us by email (isi@isi.it) or write to:
#  ISI Foundation, Via Alassio 11/c, 10126 Torino, Italy.
#
#  This program was written by André Panisson <panisson@gmail.com>
#
import unittest
from pymobility.models import contact
from pymobility.models.contact import modelB, model_het
import numpy as np
from pymobility.models.mobility import random_waypoint, random_walk,\
    truncated_levy_walk, random_direction, gauss_markov, reference_point_group,\
    tvc

class MobilityModelTestCase(unittest.TestCase):
    
    MAX_X = 100
    MAX_Y = 100
    nr_nodes = 100
    
    def setUp(self):
        np.random.seed(0xffffff)
    
    def check_margins(self, xy, max_x, max_y):
        # check if there is a node outside the margins
        b = np.where(xy[:,0]<0)[0]
        self.assertEqual(b.size, 0, "At least a node was found outside the margins")
        b = np.where(xy[:,0]>max_x)[0]
        self.assertEqual(b.size, 0, "At least a node was found outside the margins")
        b = np.where(xy[:,1]<0)[0]
        self.assertEqual(b.size, 0, "At least a node was found outside the margins")
        b = np.where(xy[:,1]>max_y)[0]
        self.assertEqual(b.size, 0, "At least a node was found outside the margins")
            
class RandomWaypointTestCase(MobilityModelTestCase):
    
    def create_model_instance(self):
        return random_waypoint(self.nr_nodes, (self.MAX_X, self.MAX_Y), velocity=(15, 20))
    
    def test_margins(self):
        model_instance = self.create_model_instance()
        for _ in range(1000):
            xy = next(model_instance)
            self.check_margins(xy, self.MAX_X, self.MAX_Y)
    
class RandomWalkTestCase(MobilityModelTestCase):
    
    def create_model_instance(self):
        return random_walk(self.nr_nodes, (self.MAX_X, self.MAX_Y))
    
    def test_margins(self):
        model_instance = self.create_model_instance()
        for _ in range(1000):
            xy = next(model_instance)
            self.check_margins(xy, self.MAX_X, self.MAX_Y)
    
class TruncatedLevyWalkTestCase(MobilityModelTestCase):
    
    def create_model_instance(self):
        return truncated_levy_walk(self.nr_nodes, (self.MAX_X, self.MAX_Y))
    
    def test_margins(self):
        model_instance = self.create_model_instance()
        for _ in range(1000):
            xy = next(model_instance)
            self.check_margins(xy, self.MAX_X, self.MAX_Y)
    
class RandomDirectionTestCase(MobilityModelTestCase):
    
    def create_model_instance(self):
        return random_direction(self.nr_nodes, (self.MAX_X, self.MAX_Y))
    
    def test_margins(self):
        model_instance = self.create_model_instance()
        for _ in range(1000):
            xy = next(model_instance)
            self.check_margins(xy, self.MAX_X, self.MAX_Y)
    
class GaussMarkovTestCase(MobilityModelTestCase):
    
    def create_model_instance(self):
        return gauss_markov(self.nr_nodes, (self.MAX_X, self.MAX_Y))
    
    def test_margins(self):
        model_instance = self.create_model_instance()
        for _ in range(1000):
            xy = next(model_instance)
            self.check_margins(xy, self.MAX_X, self.MAX_Y)
    
class ReferencePointGroupTestCase(MobilityModelTestCase):
    
    def create_model_instance(self):
        return reference_point_group(self.nr_nodes, (self.MAX_X, self.MAX_Y))
    
    def test_margins(self):
        model_instance = self.create_model_instance()
        for _ in range(1000):
            xy = next(model_instance)
            self.check_margins(xy, self.MAX_X, self.MAX_Y)
    
class TimevariantCommunityTestCase(MobilityModelTestCase):
    
    def create_model_instance(self):
        return tvc(self.nr_nodes, (self.MAX_X, self.MAX_Y))
    
    def test_margins(self):
        model_instance = self.create_model_instance()
        for _ in range(1000):
            xy = next(model_instance)
            self.check_margins(xy, self.MAX_X, self.MAX_Y)

class ContactModelTestCase(unittest.TestCase):
    
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

    def test_dynamic_gnp(self):
        n = 10
        dgnp = contact.dynamic_gnp(n, 1.)
        for i in range(5):
            contacts = next(dgnp)
            self.assertEqual(len(contacts), (n*(n-1))/2)
            for i,j in contacts:
                self.assertNotEqual(i,j)
                
    def test_dynamic_gnm(self):
        n = 10
        m = 5
        dgnm = contact.dynamic_gnm(n, m)
        for i in range(5):
            contacts = next(dgnm)
            self.assertEqual(len(contacts), m)
            for i,j in contacts:
                self.assertNotEqual(i,j)
                
    def test_edge_markovian_empty(self):
        n = 10
        p = 0.
        q = 1.
        g = 0.
        em = contact.edge_markovian(n, p, q, g)
        for i in range(5):
            contacts = next(em)
            self.assertEqual(len(contacts), 0)
            for i,j in contacts:
                self.assertNotEqual(i,j)
                
    def test_edge_markovian_alternating(self):
        n = 10
        p = 1.
        q = 1.
        g = 0.
        em = contact.edge_markovian(n, p, q, g)
        for i in range(5):
            contacts = next(em)
            if i%2 == 0:
                self.assertEqual(len(contacts), (n*(n-1))/2)
            else:
                self.assertEqual(len(contacts), 0)
            for i,j in contacts:
                self.assertNotEqual(i,j)

    def test_continuous_time_edge_markovian(self):
        n = 10
        lmbd = 50.
        model = contact.continuous_time_edge_markovian(n, lmbd)
        for i in range(50):
            contacts = next(model)
            
    def test_broad_continuous_time_edge_markovian(self):
        n = 10
        alpha = 1.6
        model = contact.broad_continuous_time_edge_markovian(n, 1.6)
        for i in range(50):
            contacts = next(model)

