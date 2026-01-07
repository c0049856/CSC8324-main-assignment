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
    
    for i in range(1, len(selections)):
        series = result[:,i]
        name = selections[i]
        p.plot(times, series, label=str(name))
        
    p.legend()
 
    p.show()


import tellurium as te

MM = '''
model Michaelis_Menten
  
  // Species initializations:
   Ru5P = 5e-7;
   Rpe = 2e-7;
   Ru5P_Rpe = 0;    
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
   RibC = 2e-7;
   Riboflavin_RibC = 0;
   FMN = 0;
   FMN_RibC = 0;
   FAD = 0;
  
    // Reactions:
  Binding_Ru5P: Ru5P + Rpe -> Ru5P_Rpe;                 (k1 * Ru5P * Rpe);
  Dissociation_Ru5P: Ru5P_Rpe -> Ru5P + Rpe;                 (k1r * Ru5P_Rpe);
  Conversion_Ru5P: Ru5P_Rpe -> Rpe + R5P;                  (k2 * Ru5P_Rpe);
  Binding_R5P: R5P + RibA -> R5P_RibA;                 (k1 * R5P * RibA);
  Dissociation_R5P: R5P_RibA -> R5P + RibA;                 (k1r * R5P_RibA);
  Conversion_R5P: R5P_RibA -> RibA + D2B4P;               (k2 * R5P_RibA);
  Binding_D2B4P: D2B4P + RibH -> D2B4P_RibH;             (k1 * D2B4P * RibH);
  Dissociation_D2B4P: D2B4P_RibH -> D2B4P + RibH;             (k1r * D2B4P_RibH);
  Conversion_D2B4P: D2B4P_RibH -> RibH + D8RL;              (k2 * D2B4P_RibH);
  Binding_D8RL: D8RL + RibE -> D8RL_RibE;               (k1 * D8RL * RibE);
  Dissociation_D8RL: D8RL_RibE -> D8RL + RibE;               (k1r * D8RL_RibE);
  Conversion_D8RL: D8RL_RibE -> RibE + Riboflavin;         (k2 * D8RL_RibE);
  Binding_Riboflavin: Riboflavin + RibC -> Riboflavin_RibC;   (k1 * Riboflavin * RibC);
  Dissociation_Riboflavin: Riboflavin_RibC -> Riboflavin + RibC;   (k1r * Riboflavin_RibC);
  Conversion_Riboflavin: Riboflavin_RibC -> RibC + FMN;          (k2 * Riboflavin_RibC);
  Binding_FMN: FMN + RibC -> FMN_RibC;                 (k1 * FMN * RibC);
  Dissociation_FMN: FMN_RibC -> FMN + RibC;                 (k1r * FMN_RibC);
  Conversion_FMN: FMN_RibC -> RibC + FAD;                 (k2 * FMN_RibC);
  Biomass: Riboflavin ->;                   (k_d * Riboflavin);
  Synthesis_Ru5P: -> Ru5P; (k_s);
  Synthesis_R5P: -> R5P; (k_sR5P);

  // Variable initialization:
  k1  = 1e6; 
  k1r = 1e-4; 
  k2  = 0.1;
  
  k_s = 0;
  k_d = 0;
  k_sR5P = 0;
end'''

r = te.loada(MM)
result = r.simulate(0, 100, 1000)

# plot the simulation output 
#r.plot()
my_plot(r, result, sizeX=12, sizeY=8)

# export to SBML
r.exportToSBML('MMforFluxALL.xml', current = False)



