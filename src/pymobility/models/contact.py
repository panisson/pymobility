# coding: utf-8
#
#  Copyright (C) 2008-2010 Istituto per l'Interscambio Scientifico I.S.I.
#  You can contact us by email (isi@isi.it) or write to:
#  ISI Foundation, Viale S. Severo 65, 10133 Torino, Italy. 
#
#  This program was written by André Panisson <panisson@gmail.com>
#
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
