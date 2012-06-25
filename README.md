Description
-----------
This is a python library for simulation of mobility and contact models.

The following mobility models that simulate node movements are available:

- Random Walk
- Random Waypoint
- Random Direction
- Truncated Levy Walk
- Gauss-Markov
- Reference Point Group
- Time-variant Community

In addition to models that simulate node position in time, this library also provides some models 
that simulate node contacts:
- The Random Contact model, that simulate random contacts in time
- The Model B, as described in [1]

Dependencies
------------
NumPy and Matplotlib

Examples
--------
To run the following examples, clone this repository, open a shell prompt, go to the src directory 
and run the python shell:
```bash
export PYTHONPATH=. && python
```
### Mobility Models
Once you are in the python shell, you can instantiate the mobility models available 
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
For example, to create a Random Contact model instance, use the following commands:
```python
>>> from pymobility.models.contact import random_contact
>>> rc = random_contact(200, avg_contacts=5)
```
This will create a Random Contact model instance with 200 nodes that generates a list of contacts
at each time step. 
The average number of contacts generated in each step is defined by the *avg_contacts* parameter.
Now you can iterate over the list of contacts produced in each time step:
```python
>>> for contacts in rc:
...     for source, target in contacts:
...         print (source, target)
```

### Simulation and Visualization
The script pymobility/simulation.py, under the src directory, contains examples on how to run different models 
and to plot the points in a simulation area using Matplotlib.
To run it, open a shell prompt, go to the src directory and run the following command:
```bash
export PYTHONPATH=. && python pymobility/simulation.py
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
[1] Stehle, J., Barrat, A. & Bianconi, G. Dynamical and bursty interactions in social networks. Physical Review E 81, 1-4 (2010). Available online at: http://link.aps.org/doi/10.1103/PhysRevE.81.035101
