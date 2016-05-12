###################################################
##  
##  Perform fits to the UCLA measurements (ELOG 77).
##  The fit parameters will be implemented in IceMC.
##
###################################################
import numpy, scipy, pylab
from scipy.optimize import curve_fit

clrs=['k','b','r','c','fuchsia','goldenrod','darkslategray','orchid','sienna','indianred','m']

def snell(inc):
  return scipy.arcsin((1./1.5)*scipy.sin(inc))

# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))



ddir = '/Users/michael/devl/icemc/trunk/data/roughness/'
df = numpy.loadtxt(ddir + 'masterfile.txt')
grit = df[:,0]
t0 = df[:,1]
t = df[:,2] # positive values are TOWARDS THE SURFACE NORMAL
p = df[:,3]


coeff_array=[]

## ELOG 77, Table 4: Measured fraction of light transmitted though flat glass at various angles of incidence compared to the theoretically predicted values. The laser power is 136 microW.
## IncidenceAngle, PowerTransmitted, FracTransmitted_meas, FresnelFactor, LossFactor, FracTransmitted_theory 
table4 = [[ 0, 114, 0.84, 0.92, 0.91, 0.84], \
          [10, 111, 0.82, 0.92, 0.91, 0.84], \
          [20, 110, 0.81, 0.92, 0.91, 0.84], \
          [30, 108, 0.79, 0.91, 0.91, 0.83], \
          [40, 103, 0.76, 0.90, 0.90, 0.81], \
          [50,  96, 0.71, 0.87, 0.90, 0.78], \
          [60,  83, 0.61, 0.80, 0.89, 0.72], \
          [70,  61, 0.45, 0.65, 0.89, 0.58]]


grit_uniq = numpy.unique(grit)

f=pylab.figure(num=0,figsize=(7,6))
for G in grit_uniq:
  f.clf()
  f.suptitle('Grit %i'%G, fontsize=20)
  t0_uniq = numpy.unique(t0[numpy.where(grit==G)])
  i=0
  for T in t0_uniq:
    q = numpy.where( (grit==G)&(t0==T) )
    # measurements
    pylab.plot(t[q], p[q], 'd-', color=clrs[i], label='$\Theta_{0}=%i$'%T)
    # gaussian fit
    p0 = [max(p[q]), scipy.mean(t[q]), scipy.std(t[q])]
    coeff, var_matrix = curve_fit(gauss, t[q], p[q], p0=p0, maxfev=1000000)
    X = numpy.linspace(min(t[q]), max(t[q]), 1000)
    yfit = gauss(X, *coeff)
    pylab.plot(X, yfit, '--', color=clrs[i])
    i += 1
    coeff_array.append( [G,T,coeff[0],coeff[1],coeff[2]] )
    print var_matrix
  pylab.legend(loc=0,numpoints=1)
  pylab.grid()
  pylab.xlabel('$\Theta$ [deg]')
  pylab.ylabel('Power [$\mu$W]')
  pylab.draw()
  pylab.savefig(ddir+'plot_powerdist_grit%i.png'%G)




coeff_array = numpy.transpose(numpy.transpose(coeff_array))

f.clf()
q1 = numpy.where(coeff_array[:,0]==400)
pylab.plot(coeff_array[:,1][q1], coeff_array[:,2][q1],'k.', label='400 grit')
q1 = numpy.where(coeff_array[:,0]==1000)
pylab.plot(coeff_array[:,1][q1], coeff_array[:,2][q1],'r.', label='1000 grit')
q1 = numpy.where(coeff_array[:,0]==1500)
pylab.plot(coeff_array[:,1][q1], coeff_array[:,2][q1],'b.', label='1500 grit')
pylab.legend(loc=0,numpoints=1)
pylab.grid()
pylab.xlabel('$\Theta_{0}$ [deg]')
pylab.ylabel('Fit Amplitude [$\mu$W]')
pylab.draw()
pylab.savefig(ddir+'plot_fit_amplitude.png')

f.clf()
q1 = numpy.where(coeff_array[:,0]==400)
pylab.plot(coeff_array[:,1][q1], coeff_array[:,3][q1],'k.', label='400 grit')
q1 = numpy.where(coeff_array[:,0]==1000)
pylab.plot(coeff_array[:,1][q1], coeff_array[:,3][q1],'r.', label='1000 grit')
q1 = numpy.where(coeff_array[:,0]==1500)
pylab.plot(coeff_array[:,1][q1], coeff_array[:,3][q1],'b.', label='1500 grit')
pylab.legend(loc=0,numpoints=1)
pylab.grid()
pylab.xlabel('$\Theta_{0}$ [deg]')
pylab.ylabel('Fit $\mu$ [deg]')
pylab.draw()
pylab.savefig(ddir+'plot_fit_mean.png')

f.clf()
q1 = numpy.where(coeff_array[:,0]==400)
pylab.plot(coeff_array[:,1][q1], coeff_array[:,4][q1],'k.', label='400 grit')
q1 = numpy.where(coeff_array[:,0]==1000)
pylab.plot(coeff_array[:,1][q1], coeff_array[:,4][q1],'r.', label='1000 grit')
q1 = numpy.where(coeff_array[:,0]==1500)
pylab.plot(coeff_array[:,1][q1], coeff_array[:,4][q1],'b.', label='1500 grit')
pylab.legend(loc=0,numpoints=1)
pylab.grid()
pylab.xlabel('$\Theta_{0}$ [deg]')
pylab.ylabel('Fit $\sigma$ [deg]')
pylab.draw()
pylab.savefig(ddir+'plot_fit_sigma.png')