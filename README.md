Description
-----------
This is a python library to simulation of mobility and contact models.

The following mobility models that simulate node movements are available:

- Random Walk
- Random Waypoint
- Random Direction
- Truncated Levy Walk
- Gauss-Markov
- Reference Point Group
- Time-variant Community

In addition to models that simulate node position in time, this library also provides some models that simulate node contacts:
- The Random Contact model, that simulate random contacts in time
- The Model B, as described in [1]

Dependencies
------------
NumPy and Matplotlib

Examples
--------
The pymobility/simulation.py script, under the src directory, contains some examples on how to run each model.
To run it on a Linux shell, install the dependencies, go to the src directory and run the following command:
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
