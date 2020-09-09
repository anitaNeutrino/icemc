### Quick program to plot already produced spectra from different scenarios

import numpy as np
import pylab as pl
import os
# infinite print outls -lrth
np.set_printoptions(threshold=np.inf)
# Scale font
pl.rcParams.update({'font.size': 22})

envName='ICEMC_SRC_DIR'
subDir='/fluxes/'
if envName in os.environ:
    sourceDir = os.environ.get(envName)

number = 1000000
sourceEvolution = 1
cutoff = 1

fig = pl.figure(1)
plE = fig.add_subplot(111)
plE.margins(x=0)
#plE.set_ylim([1.3*10**-8, 6*10**-8])
plE.set_xlim([10**9, 10**12])
#plE.set_ylim([10**-10.5, 10**-5])
#plE.set_xlim([10**5, 10**11])
plE.set_title('Neutrino flux spectrum')
#plE.set_xlabel(r'$\log_{10}(E/{\rm eV})$')
plE.set_xlabel(r'$E{\rm (GeV)}$')
plE.set_ylabel(r'$E^{2} J (\rm GeV cm^{-2} sr^{-1} s^{-1})$')
plE.semilogx()
plE.semilogy()

#emaxRange = [20.1, 21.1, 23.1]
emaxRange = [20.0, 21.0, 23.0]
#alphaRange = [2.0, 2.5, 2.9]
#alphaRange = [1.0]
alpha = 2.0
cutoffRange = [1]
szRange = [1.5]
#szRange = [1.5,3.0]

for sourceEvolutionRedshiftCap in szRange:
    if(sourceEvolutionRedshiftCap == szRange[0]):
        colorType='Crimson'
        style=''
    #if(sourceEvolutionRedshiftCap == szRange[1]):
        #style=':'
        #colorType='Turquoise'
    #if(alpha == alphaRange[2]):
        #colorType='Blue'
    for emax in emaxRange:
        data = np.loadtxt(sourceDir + subDir + 'simFlux_E' + str(emax) + '_A' + str(alpha) + '_N' + str(number) + '_c' + str(cutoff) + '_s' + str(sourceEvolution) + '_z' + str(sourceEvolutionRedshiftCap) + '.dat', delimiter='  ', dtype=float, skiprows=1, usecols=range(2))
        #data = np.loadtxt(sourceDir + subDir + 'simFlux_E' + str(emax) + '_A' + str(alpha) + '_N' + str(number) + '_c' + str(cutoff) + '_s' + str(sourceEvolution) + '.dat', delimiter='  ', dtype=float, skiprows=1, usecols=range(2))
        #data = np.loadtxt(sourceDir + subDir + 'simFlux_E' + str(emax) + '_A' + str(alpha) + '_N' + str(number) + '_c' + str(cutoff) + '_s' + str(sourceEvolution) + '.dat', delimiter='  ', dtype=float, skiprows=1, usecols=range(2))
        if(emax == emaxRange[0]):
            style=''
        if(emax == emaxRange[1]):
            style='--'
        if(emax == emaxRange[2]):
            style=':'
        # convert to different units for comparison
        data[:,0]=10**(data[:,0]-9) # logeV -> GeV
        data[:,1]=10**(data[:,1]) # log-> normal
        #plE.plot(data[:,0], data[:,1], style, color=colorType, lineWidth=4.0, label='Emax = 10' + '$^{'+str(emax)+'}$' + 'eV,' + r' $\alpha$ =' + str(alpha))
        #plE.plot(data[:,0], data[:,1], style, color=colorType, lineWidth=4.0, label='sourceEvolutionRedshiftCap ' + str(sourceEvolutionRedshiftCap))
        plE.plot(data[:,0], data[:,1], style, color=colorType, lineWidth=4.0, label='Emax = 10' + '$^{'+str(emax)+'}$' + 'eV')
pl.legend()
pl.show()
