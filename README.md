# MReaDy (Multiprocess Reaction Dynamics program) 

MReaDy is a program for studying reactive dynamic systems using a global potential energy surface (gPES). 
Potential Energy Surfaces (PES), based on ab initio calculations,
is a powerful tool to study the rate of elementary reactions and their dynamics,
being useful to compute state-to-state rate constants.
In a more complex mechanism, we will be in the presence of different and simultaneous elementary reactions,
corresponding to all the possible reactive and non-reactive collisions between the species present and leading to the respective products.
Attempting to build a traditional PES for such a system easily becomes impossible.
To circumvent this problem, a global Potential Energy Surface (gPES) can be defined by integrating various PESs,
each representing an elementary reaction expected to play a role in the chemical process. 
MReaDy is built in such a way and performs reactive dynamic calculations based on such gPES.
There are several versions of MReaDy, and gradually they will be uploaded.
We present MReaDy Base version since most of its working principles are common to the remaining versions.
It simulates a hydrogen combustion with molecules that go up to 4 atoms. 
Working principles are presented on the reference given.

How to reference MReaDy:
César Mogo and João Brandão
J. Comput. Chem. 2014, 35, 1330-1337.
DOI: 10.1002/jcc.23621
