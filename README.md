Description
-----------
This is a python library for simulation of mobility and contact models.

The following mobility models that simulate node movements are available:

- Random Walk
- Random Waypoint
- Random Direction
- Truncated Levy Walk [[1]](#references)
- Gauss-Markov [[2]](#references)
- Reference Point Group Mobility model [[3]](#references)
- Time-variant Community [[4]](#references)

In addition to models that simulate node position in time, this library also provides some models 
that simulate node contacts:
- Dynamic G(n,p), a dynamic version of the Erdős–Rényi model [[5]](#references)
- Dynamic G(n,m), another dynamic version of the Erdős–Rényi model [[5]](#references)
- The edge-Markovian dynamic graph model [[6]](#references)
- The continuous-time edge-Markovian dynamic graph model [[7]](#references)
- The Model B, as described in [[8]](#references)

Installation
------------

setuptools - from Git repository

```bash
> git clone git://github.com/panisson/pymobility.git
> cd pymobility
> python setup.py install (run as admin/root)
```

Dependencies
------------
NumPy and Matplotlib

Examples
--------
### Mobility Models
In a python shell, you can instantiate the mobility models available 
in the pymobility.models.mobility package.
For example, to create a Random Waypoint model instace, use the following commands:
```python
>>> from pymobility.models.mobility import random_waypoint
>>> rw = random_waypoint(200, dimensions=(100, 100), velocity=(0.1, 1.0), wt_max=1.0)
```
This will create a Random Waypoint instance with 200 nodes in a simulation area of 100x100 units, 
velocity chosen from a uniform distribution between 0.1 and 1.0 units/step
and maximum waiting time of 1.0 steps.
This object is a generator that yields the position of the nodes in each step.
For example, to print a 2-dimensional array with the node positions produced in the first step, call
```python
>>> positions = next(rw)
>>> print positions
```
You can also iterate over the positions produced in each step:
```python
>>> for positions in rw:
...     print positions
... 
```
### Contact Models
The contact models are available in the pymobility.models.contact package.
Differently from the mobility package that yield node positions for each time step, 
the contact models yield a list of node pairs that were in contact at each time step.
For example, to create a Dynamic G(n,p) model instance, use the following commands:
```python
>>> from pymobility.models.contact import dynamic_gnp
>>> m = dynamic_gnp(200, p=0.01)
```
This will create a Dynamic G(n,p) model instance with 200 nodes that generates a list of contacts
at each time step. 
The probability of any edge appearing at each step is defined by the parameter *p*.  
Now you can iterate over the list of contacts produced in each time step:
```python
>>> for contacts in m:
...     for source, target in contacts:
...         print (source, target)
```

### Simulation and Visualization
The script pymobility/simulation.py, under the src directory, contains examples on how to run different models 
and to plot the points in a simulation area using Matplotlib.
To run it, open a shell prompt, go to the src directory and run the following command:
```bash
python pymobility/simulation.py
```

Contributing
------------
If you have a Github account please fork the repository,
create a topic branch, and commit your changes.
Then submit a pull request from that branch.

License
-------
Written by André Panisson <panisson@gmail.com>  
Copyright (C) 2012 Istituto per l'Interscambio Scientifico I.S.I.  
You can contact us by email (isi@isi.it) or write to:  
ISI Foundation, Via Alassio 11/c, 10126 Torino, Italy.  

The Model B was implemented with the contribution of Juliette Stehle.

References
----------
[1] Injong Rhee, Minsu Shin, Seongik Hong, Kyunghan Lee, and Song Chong. On the Levy-Walk Nature of Human Mobility. 
    In 2008 IEEE INFOCOM - Proceedings of the 27th Conference on Computer Communications, pages 924-932. April 2008.

[2] Camp, T., Boleng, J. & Davies, V. A survey of mobility models for ad hoc network research. 
    Wireless Communications and Mobile Computing 2, 483-502 (2002).

[3] Xiaoyan Hong, Mario Gerla, Guangyu Pei, and Ching-Chuan Chiang. 1999. 
    A group mobility model for ad hoc wireless networks. In Proceedings of the 
    2nd ACM international workshop on Modeling, analysis and simulation of 
    wireless and mobile systems (MSWiM '99). ACM, New York, NY, USA, 53-60.

[4] Wei-jen Hsu, Thrasyvoulos Spyropoulos, Konstantinos Psounis, and Ahmed Helmy, 
    Modeling Time-variant User Mobility in Wireless Mobile Networks, INFOCOM 2007, May 2007.

[5] Andrea E. F. Clementi, Francesco Pasquale, Angelo Monti, and Riccardo Silvestri. 2007. 
    Communication in dynamic radio networks. In Proceedings of the twenty-sixth annual 
    ACM symposium on Principles of distributed computing (PODC '07). ACM, New York, NY, USA, 205-214.

[6] Andrea E.F. Clementi, Claudio Macci, Angelo Monti, Francesco Pasquale, and Riccardo Silvestri. 2008. 
    Flooding time in edge-Markovian dynamic graphs. In Proceedings of the 
    twenty-seventh ACM symposium on Principles of distributed computing (PODC '08). 
    ACM, New York, NY, USA, 213-222.

[7] Augustin Chaintreau, Abderrahmen Mtibaa, Laurent Massoulie, and Christophe Diot. 2007. 
    The diameter of opportunistic mobile networks. In Proceedings of the 
    2007 ACM CoNEXT conference (CoNEXT '07). ACM, New York, NY, USA, , Article 12 , 12 pages.

[8] Stehle, J., Barrat, A. & Bianconi, G. Dynamical and bursty interactions in social networks. 
    Physical Review E 81, 1-4 (2010). Available online at: http://link.aps.org/doi/10.1103/PhysRevE.81.035101
