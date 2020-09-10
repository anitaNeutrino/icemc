#! /usr/bin/env python3
############################################################################
# Desc:
# 1D simulation for producing neutrino flux spectra at ultra-high energies
# Interactions include photopion production and nuclear decay
# Both nucleon and neutrino fluxes are produced, normalized, and plotted
# Spectra are then output into the standard icemc input spectrum format
############################################################################
# Basic usage:
# python simulatedNeutrinoFlux.py
#
# For options and explanations: 
# python simulatedNeutrinoFlux.py -h
#
# Example with some paramaters
# python simulatedNeutrinoFlux.py -e 23.1 a = 2.2 -n 300000 -d 1
# This sets maximum produced energy to 10^23.1 EeV, alpha (power law index) 
# to 2.2, the number of CRs propagated to 300000, and draws neutrino spectra
############################################################################


### Libs, crpropa
from crpropa import *

import numpy as np
from numpy import inf
import pylab as pl
import os
import argparse
import sys

### Settings
# infinite print out
np.set_printoptions(threshold=np.inf)
# Scale font
pl.rcParams.update({'font.size': 22})

### Input options
print ('_____INPUT OPTIONS_____')
parser = argparse.ArgumentParser('InputParams')
parser.add_argument('-e', '--emax', default=21, help='desc: Maximum energy (eV) exponent of cosmic rays from the source <float> Default = 23.1 (10^23.1) <conditions> High energy (~ZeV scale)')
parser.add_argument('-a', '--alpha', default=2.0, help='desc: Source power law index value <float> Default = 2.0')
parser.add_argument('-n', '--number', default=1000000, help='desc: Number of cosmic rays to propagate <int> Default = 1000000')
parser.add_argument('-d', '--draw', default=1, help='desc: Draw spectra <int> Default = 1 <conditions> Off = 0, On = 1, Show all spectra = 2')
parser.add_argument('-c', '--cutOffOption', default=0., help='desc: Add exponential cuf-off factor to injection spectrum <int> Default = 0 <conditions> Off = 0, On = 1')
parser.add_argument('-s', '--sourceEvolutionParam', default=3, help='desc: Source evolution param (sometimes called m) <int> Default = 3 <float> Off = 0, otherwise use that exponent')
parser.add_argument('-S', '--scratchDir', default='.', help='desc: Scratch directory (for writing temporary outputs) <int> Default = .)')
parser.add_argument('-k', '--keep', default=False, help='desc Keep text output files (these can be massive and will likely get overwritten anyway) <bool> Default = False)')
parser.add_argument('-b', '--batchMode', default=0, help='desc: Enable batch mode to run on clusters <int> Default = 1 <conditions> Off = 0, On = 1')
parser.add_argument('-z', '--sourceEvolutionRedshiftCap', default=1.5, help='desc: Source evolution scales with (1+z)^m until this redshift value. NOTE: This is the not maximum redshift, just until source evolution plateaus! <float> Default = 1.5')
parser.add_argument('-m', '--minRedShift', default=0, help='desc: the minimum redshift! <float> Default = 0')
parser.add_argument('-M', '--maxRedShift', default=4, help='desc: the maximum redshift! <float> Default = 4')
parser.add_argument('-o', '--output', default='', help='desc: output file name (otherwise auto-generated from optionss)')
parser.add_argument('-D', '--outputDir', default='', help='desc: output dir (otherwise ICEMC_SRC_DIR or cwd)')
args = parser.parse_args()
emax = float(args.emax)
alpha = float(args.alpha)
number = int(args.number)
draw = int(args.draw)
cutOffOption = int(args.cutOffOption)
sourceEvolution = int(args.sourceEvolutionParam)!=0
batchMode = int(args.batchMode)
sourceEvolutionRedshiftCap = float(args.sourceEvolutionRedshiftCap)

print('Maximum energy = 10^' + str(emax) + ' eV')
if(emax < 20):
    sys.exit('Specify a higher MAXIMUM possible energy for cosmic rays from the source. Typical values 10^20 - 10^23 eV')
print('Alpha = ' + str(alpha))
print('Source evolution redshift cap = ' + str(sourceEvolutionRedshiftCap))
print('Cosmic rays simulated from the source = ' + str(number))
if(draw == 1 or draw == 2):
    print('Plotting spectra')
if(cutOffOption==1):
    print('Applying cut-off to injection spectrum')
if(sourceEvolution):
    print('Applying source evolution')
if(batchMode==1):
    print('Running in batch mode')
    # Batch farm output
    pl.switch_backend('agg') # required output for batch farm (none)
else:
    print('Running in normal mode')
print ('_______________________')
### End input options

######## Start CRPROPA assumptions and settings
# Only interested in neutrinos as secondaries
neutrinos = True
photons = False
electrons = False

# Set up variables
### NOTE, these are the inner and outer BIN EDGES, respectively
### The actual energies used are the bin CENTERS, i.e. 16.95, 21.05
energyMin = 10.**16.9
energyMax = 10.**emax # range: 100 EeV -> 10^5 EeV
logEmin = np.log10(energyMin) # log of min energy
logEmax = np.log10(energyMax) # log of max energy
redshiftMin = float(args.minRedShift);
redshiftMax = float(args.maxRedShift); 
mParam = float(args.sourceEvolutionParam); 

# Set up modules
m = ModuleList()
m.add(SimplePropagation(10*kpc, 10*Mpc))
m.add(Redshift())
m.add(PhotoPionProduction(CMB, photons, neutrinos))
m.add(PhotoPionProduction(IRB, photons, neutrinos))
m.add(NuclearDecay(electrons, photons, neutrinos))
m.add(ElectronPairProduction(CMB))
m.add(ElectronPairProduction(IRB))
m.add(MinimumEnergy(energyMin * eV))

nucleonOutputFile = args.scratchDir + '/out_nucleons_E' + str(emax) + '_A' + str(alpha) + '_N' + str(number) + '_c' + str(cutOffOption) + '_s' + str(sourceEvolution) + '_z' + str(sourceEvolutionRedshiftCap) + '.txt'
neutrinoOutputFile = args.scratchDir + '/out_neutrinos_E' + str(emax) + '_A' + str(alpha) + '_N' + str(number) + '_c' + str(cutOffOption) + '_s' + str(sourceEvolution) + '_z' + str(sourceEvolutionRedshiftCap) + '.txt'

# Observer for nucleons
obs1 = Observer()
obs1.add(ObserverPoint())
obs1.add(ObserverNeutrinoVeto())  # we don't want neutrinos here
output1 = TextOutput(nucleonOutputFile, Output.Event1D)
obs1.onDetection(output1)
m.add(obs1)

# Observer for neutrinos
obs2 = Observer()
obs2.add(ObserverPoint())
obs2.add(ObserverNucleusVeto())  # we don't want nucleons here
output2 = TextOutput(neutrinoOutputFile, Output.Event1D)
obs2.onDetection(output2)
m.add(obs2)

# Source: protons with power-law spectrum from uniformly distributed sources with redshift z = 0-3
source = Source()
#if(sourceEvolution==0):
source.add(SourceUniform1D(redshift2ComovingDistance(redshiftMin), redshift2ComovingDistance(redshiftMax)))
source.add(SourceRedshift1D())

if sourceEvolution:
    source.add(SourceRedshiftEvolution(mParam, 0, sourceEvolutionRedshiftCap)) #Source evolution
# E^alpha spec (can't handle exponential cut off, so this is accounted for later)
source.add(SourcePowerLawSpectrum(energyMin * eV, energyMax * eV, -1 * alpha))
#source.add(SourcePosition(100 * Mpc))
source.add(SourceParticleType(nucleusId(1, 1))) # Protons only for now!

# Run simulation for primaries and propagate all secondaries
print('Running CRPropa...')
m.setShowProgress(True)
m.run(source, number, True)
print('Finished running CRPropa')

output1.close()
output2.close()

# Energies in EeV
nucleons = np.genfromtxt(nucleonOutputFile, names=True) 
neutrinos = np.genfromtxt(neutrinoOutputFile, names=True)

# Don't need input files from here on, we re-generate them every time program is used
if not args.keep: 
  os.remove(nucleonOutputFile)
  os.remove(neutrinoOutputFile)

####
print('Normalizing, generating flux spectra...')
####### Overall figure
fig = pl.figure(1)
countNorm = float(len(nucleons['E']))
countNormNeutrino = float(len(neutrinos['E']))
#### Nucleons
### Hist setup
energyStep = 0.10 # how spaced our log
lE = pl.log10(nucleons['E']) + 18.  # energy in log10(E/eV))
lEbins = pl.arange(logEmin, logEmax, energyStep)  # logarithmic bins
lEcens = (lEbins[1:] + lEbins[:-1]) / 2  # logarithmic bin centers

## We apply a posteriori weight binning when filling the histograms
## as have the original energy E0 available from the output files
if (cutOffOption==1):
    exponentialCutoffN = np.exp(- 1 * (nucleons['E0'] * 10**18)/(energyMax))
    exponentialCutoffV = np.exp(- 1 * (neutrinos['E0'] * 10**18)/(energyMax))
else:
    exponentialCutoffN = nucleons['E0']/nucleons['E0'] # just so we have array of 1 same length as data for weights
    exponentialCutoffV = neutrinos['E0']/neutrinos['E0'] # see above
dE = 10**lEbins[1:] - 10**lEbins[:-1]  # bin widths
JEE = pl.histogram(lE, weights=exponentialCutoffN, bins=lEbins)[0] / dE * 10**lEcens * 10**lEcens

### Nucleon number hist
nNuc = fig.add_subplot(241)
nNuc.margins(x=0)
    
nNuc.hist(lE, weights=exponentialCutoffN, bins=lEbins, histtype='step', log=True, density=False, color='Blue')
nNuc.set_title('Nucleon count vs energy')
nNuc.set_xlabel(r'$\log_{10}(E/{\rm eV})$')
nNuc.set_ylabel(r'$n_{\rm nucleon}$ [a.u.]')

### Nucleon JEE hist
jEENuc = fig.add_subplot(242)
jEENuc.margins(x=0)
#pl.plot(lEcens, J,  color='SaddleBrown')
jEENuc.plot(lEcens, JEE, color='Blue')
jEENuc.set_title('Nucleon flux spectrum')
jEENuc.set_xlabel(r'$\log_{10}(E/{\rm eV})$')
jEENuc.set_ylabel(r'$E^{2} J$ [a.u.]')
jEENuc.semilogy()
jEENuc.set_xlim([17, 21])

########## AUGER
J0 = 3.30 * 10**-19.  #eV^-1 km^-2 sr^-1 yr^-1
### Conversions for plots
yearToSecond = 3.15*10**7.  #unit conversion of y axis in fig3 in auger paper (earlier, if needed for comparison)
kmToM = 1000.
kmToCM = 100000.
eVToGeV = 10**-9
#J0*=(yearToSecond)**(-1) * (kmToM)**(-2)

#E < Eankle
EAnkle = 4.82 *10**18. #eV
gamma1 = 3.29

#E >= Eankle
Es = 42.09 *10**18. #eV
gamma2 = 2.60
deltaGamma = 3.14

## Auger function
def JEEA(E):   #Below ankle                                                                                #Above ankle
    return (E<EAnkle)*(J0*(E/EAnkle)**(-gamma1) * E * E) + (E>=EAnkle)*(J0*(E/EAnkle)**(-gamma2) * (1 + (EAnkle/Es)**deltaGamma) * (1 + (E/Es)**deltaGamma)**(-1) * E * E)
#E1 = np.arange(energyMin, energyMax, 10**16)
## Align to bin width
#logE1 = np.arange(logEmin+0.05, logEmax-0.05, 0.10)
#print(logE1[21]) # <- bin associated with 10**19

#print(logE1)
#print(lEcens)

#Auger
#augerSpec = fig.add_subplot(243)
#augerSpec.margins(x=0)
#augerSpec.plot(lEcens, JEEA(10**logE1), 'k') ## lEcens = logE1
#augerSpec.set_yscale('log')
#augerSpec.set_title('AUGER flux spectrum vs energy')
#augerSpec.set_xlabel(r'$\log_{10}(E/{\rm eV})$')
#augerSpec.set_ylabel(r'$E^{2} J (\rm eV km^{-2} sr^{-1} yr^{-1})$')

##################### GET AUGER DATA FLUX AT A SPECIFIC ENERGY: NORMALIZATION
###### At 10^logEnergyNorm eV, fluxNorm = flux(E = logEnergyNorm).
###### Thus we scale the bins, such that
###### Find out which bin relates to energy 10^19.55 (changes due to eMin)
energyList = np.array(lEcens).ravel().tolist()
energyListR = list(np.around(np.array(energyList),2))
logEnergyNorm= 19.55 # Often used in papers
energyNorm = 10.**logEnergyNorm # Energy to normalize to
JEEAtNorm = JEEA(energyNorm)

normIteration=0
normIndex = energyListR.index(logEnergyNorm) # <- index associated with 10**energyNorm when starting from our minimum
fixedBinJEE = JEE[normIndex]

while(fixedBinJEE == 0):
    if(normIteration==0):
        print('Normalization not possible at 10^' + str(logEnergyNorm) + 'eV. This is because the nucleon spectrum has no counts in this bin! Normalizing at a lower energy instead... It is highly recommended that you run again with higher statistics.')
    normIteration+=1
    logEnergyNorm-= energyStep # our minimum energy step
    energyNorm = 10.**logEnergyNorm # Energy to normalize to
    normIndex = energyListR.index(np.around(logEnergyNorm,2)) # <- index associated with 10**energyNorm when starting from our minimum
    fixedBinJEE = JEE[normIndex]
    JEEAtNorm = JEEA(energyNorm)
   
print('Normalizing with respect to the AUGER spectrum at E = 10^' + str(logEnergyNorm) + 'eV')
JEEScaled = JEE * JEEAtNorm/fixedBinJEE
JEScaled = JEEScaled/(10**lEcens) # this should work

### Nucleon JEE hist
jEENucSc = fig.add_subplot(243)
jEENucSc.margins(x=0)
jEENucSc.plot(lEcens, JEEScaled, color='Blue')  #lEcens = logE1
jEENucSc.set_yscale('log')
jEENucSc.set_title('Nucleon spectrum normalized to AUGER')
jEENucSc.set_xlabel(r'$\log_{10}(E/{\rm eV})$')
jEENucSc.set_ylabel(r'$E^{2} J (\rm eV km^{-2} sr^{-1} yr^{-1})$')
jEENucSc.set_xlim([17, 21])

### Nucleon JEE hist both
jEENucScBoth = fig.add_subplot(244)
jEENucScBoth.margins(x=0)
#jEENucScBoth.plot(lEcens, JEEScaled,  c='b')
jEENucScBoth.set_title('Nucleon flux spectrum')
jEENucScBoth.set_xlabel(r'$\log_{10}(E/{\rm eV})$')
jEENucScBoth.set_ylabel(r'$E^{2} J (\rm eV km^{-2} sr^{-1} yr^{-1})$')
jEENucScBoth.semilogy()
### Want jEENuc and augerSpec on same plot
jEENucScBoth.plot(lEcens, JEEScaled, label='Simulated nucleon spectrum', color='Blue')
jEENucScBoth.plot(lEcens, JEEA(10.**lEcens), label='Auger nucleon spectrum', color='Red')
jEENucScBoth.legend(loc='upper right')
jEENucScBoth.set_xlim([17, 21])
jEENucScBoth.set_ylim([1e13, None])
############################ END NUCLEONS

############################ START NEUTRINOS

#countNormFactor = countNormNeutrino/countNorm
countNormFactor = 1

### Setup
lEv = pl.log10(neutrinos['E']) + 18  # energy in log10(E/eV))
# Binning 
JEEv = pl.histogram(lEv,weights=exponentialCutoffV, bins=lEbins)[0] / dE * 10**lEcens * 10**lEcens
JEv = pl.histogram(lEv,weights=exponentialCutoffV, bins=lEbins)[0] / dE * 10**lEcens # different scaling

### Neutrino number hist
nNeu = fig.add_subplot(245)
nNeu.margins(x=0)
nNeu.hist(lEv, weights=exponentialCutoffV, bins=lEbins, histtype='step', log=True, density=False, color='Green')
nNeu.set_title('Neutrino count vs energy')
nNeu.set_xlabel(r'$\log_{10}(E/{\rm eV})$')
nNeu.set_ylabel(r'$n_{\rm neutrino}$ [a.u.]')

#Neutrino JEE hist
jEENeu = fig.add_subplot(246)
jEENeu.margins(x=0)
jEENeu.plot(lEcens, JEEv, color='Green')
jEENeu.set_title('Neutrino flux spectrum')
jEENeu.set_xlabel(r'$\log_{10}(E/{\rm eV})$')
jEENeu.set_ylabel(r'$E^{2} J$ [a.u.]')
jEENeu.semilogy()
jEENeu.set_xlim([17, 21])

#################################

JEEvScaled = JEEv * JEEAtNorm/fixedBinJEE
JEvScaled = JEv * JEEAtNorm/fixedBinJEE

### Neutrino JEE hist scaled
jEENeuSc = fig.add_subplot(247)
jEENeuSc.margins(x=0)
jEENeuSc.set_title('Neutrino spectrum normalized to AUGER')
jEENeuSc.set_xlabel(r'$\log_{10}(E/{\rm eV})$')
jEENeuSc.set_ylabel(r'$E^{2} J (\rm eV km^{-2} sr^{-1} yr^{-1})$')
jEENeuSc.semilogy()
### Want jEENeu as final output
jEENeuSc.plot(lEcens, JEEvScaled, color='Green')
jEENeuSc.set_xlim([17, 21])

### JEE hist all
jEENeuScAll = fig.add_subplot(248)
jEENeuScAll.margins(x=0)
jEENeuScAll.set_title('Neutrino flux spectrum')
jEENeuScAll.set_xlabel(r'$\log_{10}(E/{\rm eV})$')
jEENeuScAll.set_ylabel(r'$E^{2} J (\rm eV km^{-2} sr^{-1} yr^{-1})$')
jEENeuScAll.semilogy()
### Want jEENeu and augerSpec on same plot
jEENeuScAll.plot(lEcens, JEEScaled, label='Simulated nucleon spectrum', color='Blue')
jEENeuScAll.plot(lEcens, JEEvScaled, label='Simulated neutrino spectrum', color='Green')
jEENeuScAll.plot(lEcens, JEEA(10**lEcens), label='Auger nucleon spectrum', color='Red')
jEENeuScAll.legend(loc='upper right')
jEENeuScAll.set_xlim([17, 21])
jEENeuScAll.set_ylim([1e13, None])
if(draw==1):
    pl.close()

fig2 = pl.figure(2)
########### Convert output to ice mc flux spectra
### log(E [eV])      |         E^2 J: log(E^2 J [GeV cm^-2 s^-1 sr^-1])
#Our format for E^2 J: eV km^-2 sr^-1 yr^-1
### Convert our flux 
JEEvScaledFormat = JEEvScaled * eVToGeV * kmToCM**-2 * yearToSecond**-1
JEvScaledFormat = JEvScaled * kmToCM**-2 * yearToSecond**-1
jEENeuScF = fig2.add_subplot(121)
jEENeuScF.margins(x=0)
jEENeuScF.set_title('Neutrino flux spectrum')
jEENeuScF.set_xlabel(r'$\log_{10}(E/{\rm eV})$')
jEENeuScF.set_ylabel(r'$E^{2} J (\rm GeV cm^{-2} sr^{-1} s^{-1})$')
jEENeuScF.semilogy()
jEENeuScF.plot(lEcens, JEEvScaledFormat, color='Green')
jEENeuScF.set_xlim([17, 21])

jENeuScF = fig2.add_subplot(122)
jENeuScF.margins(x=0)
jENeuScF.set_title('Neutrino flux spectrum')
jENeuScF.set_xlabel(r'$\log_{10}(E/{\rm eV})$')
jENeuScF.set_ylabel(r'$E J (\rm cm^{-2} sr^{-1} s^{-1})$')
jENeuScF.semilogy()
jENeuScF.plot(lEcens, JEvScaledFormat, color='Green')
jENeuScF.set_xlim([17, 21])

fig2 = pl.figure(2)

# If trailing have no counts in bins, cut the spectrum off there (i.e., sim ran with low data), also accounts for log(0)
JEEvSList = np.array(JEEvScaled).ravel().tolist()
if JEEvSList.count(0) > 0: # if occurences of zeroes, cut spectrum
    normIndex = JEEvSList.index(0) # find where first zero is!
    JEEvScaledFormat = JEEvScaledFormat[0:normIndex] # cut fluxes
    lEcens = lEcens[0:normIndex] # cut energies
    
    
format=np.c_[10**lEcens, JEEvScaledFormat]
fileName = args.output; 
if fileName == '':
  fileName='simFlux_E' + str(emax) + '_A' + str(alpha) + '_N' + str(number) + '_c' + str(cutOffOption) + '_s' + str(mParam) + '_z' + str(sourceEvolutionRedshiftCap) + '_m' + str(redshiftMin) + '_M' + str(redshiftMax) + '.dat'

decimalPrecision='%g'

outdir = args.outputDir

if outdir == '':
    outdir ='.' 


print('Saving spectra to ' + outdir + "/" + fileName) 

np.savetxt(outdir + '/' + fileName, format, delimiter=' ', fmt=decimalPrecision)

    # Show all plots
if(draw==1):
    print('Sim done, displaying output.')
    pl.show()

if(draw==2):
    print('Sim done, displaying output.')
    pl.show()
    
print('All done!')
