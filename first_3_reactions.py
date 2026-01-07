# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 11:48:59 2019

@author: Paolo
"""

# plotting utility from Tellurium documentation
def my_plot(r, result, sizeX=6, sizeY=4):

    import pylab as p

    if result is None:
        raise Exception("no simulation result")

    # assume result is a standard numpy array

    selections = r.timeCourseSelections

    if len(result.shape) != 2 or result.shape[1] != len(selections):
        raise Exception("simulation result columns not equal to number of selections,"
                        "likely a simulation has not been run")

    times = result[:,0]
        
    p.figure(figsize=(sizeX, sizeY))
     
    import numpy as np
    colors = p.cm.tab20(np.linspace(0, 1, len(selections)))
     
    for i in range(1, len(selections)):
         series = result[:,i]
         name = selections[i].strip('[]')
        
         key_metabolites = ["Riboflavin"]

         if str(name) in key_metabolites:
             p.plot(times, series, label=str(name),
            color=colors[i], linewidth=2.5, alpha=1.0)
         else:
             p.plot(times, series, label=str(name),
             color=colors[i], linewidth=1.0, alpha=0.6)
        
    p.legend()
    xlabel = "Time, s"; 
    ylabel = "Concentration"; 
    p.xlabel(xlabel)
    p.ylabel(ylabel)
 
    p.show()


import tellurium as te

MM = '''
model Michaelis_Menten
  
  // Species initializations:
  R5P = 5e-7;
  RibA = 2e-7;
  R5P_RibA = 0;
  D2B4P = 5e-7;
  RibH = 2e-7;
  D2B4P_RibH = 0;
  D8RL = 5e-7;
  RibE = 2e-7;
  D8RL_RibE = 0;
  Riboflavin = 0;
 
  
  // Reactions:
R5P + RibA -> R5P_RibA;         (k1 * R5P * RibA);
R5P_RibA -> R5P + RibA;         (k1r * R5P_RibA);
R5P_RibA -> RibA + D2B4P;        (k2 * R5P_RibA);
D2B4P + RibH -> D2B4P_RibH;     (k1 * D2B4P * RibH);
D2B4P_RibH -> D2B4P + RibH;     (k1r * D2B4P_RibH);
D2B4P_RibH -> RibH + D8RL;      (k2 * D2B4P_RibH);
D8RL + RibE -> D8RL_RibE;       (k1 * D8RL * RibE);
D8RL_RibE -> D8RL + RibE;       (k1r * D8RL_RibE);
D8RL_RibE -> RibE + Riboflavin; (k2 * D8RL_RibE);

  // Variable initialization:
  k1  = 1e6; 
  k1r = 1e-4; 
  k2  = 0.1;
end'''

r = te.loada(MM)
result = r.simulate(0, 175, 2000)

# plot the simulation output 
#r.plot()
my_plot(r, result, sizeX=12, sizeY=8)

# export to SBML
#r.exportToSBML('MMant.xml', current = False)



