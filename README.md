# MReaDy (Multiprocess Reaction Dynamics program) 

MReaDy is a program for studying reactive dynamic systems using a global potential energy surface (gPES). 
Potential Energy Surfaces (PES), based on ab initio calculations,
is a powerful tool to study the rate of elementary reactions and their dynamics. It is indicated to compute state-to-state rate constants.
In a more complex mechanism, we will be in the presence of different and simultaneous elementary reactions,
corresponding to all the possible reactive and non-reactive collisions between the species present and leading to the respective products.
Attempting to build a traditional PES for such a system quickly becomes impossible.
To circumvent this problem, a global Potential Energy Surface (gPES) can be defined by integrating various PESs,
each representing an elementary reaction expected to play a role in the chemical process. 
MReaDy is built in such a way that it performs reactive dynamic calculations based on such gPES.
I've explained the working principles in the reference I've given.


There are several versions of MReaDy, and they will gradually be uploaded.
A manual, file structure, and executable of the MReaDy Base version and extra files are presented here. 
It simulates a hydrogen combustion with molecules that go up to 4 atoms. 

A proper version is being prepared to submit to Github.
But if you are interested, don't hesitate to contact me at cfmogo@ualg.pt

Contributors: César Mogo, João Brandão, Carolina Rio, Wenli Wang and Daniela Coelho  


How to reference MReaDy:

César Mogo, João Brandão and D. Coelho
J. Comput. Chem. 2014, 35, 1330-1337.
DOI: 10.1002/jcc.23621
