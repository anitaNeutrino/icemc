#!/bin/python

# Modified code by Ian Best to add cross-pol data to the output as well as
# handle the second channel (With goal to provide all the data requried by
# IceMC

#Translate XF output into readable AraSim input. Jorge Torres, Nov 2017.
#You need to have all your XF output files in the same folder, with the standard
#name dipole_freq_MHz.uan



import math
import itertools
import sys
import shutil
import numpy
import glob

freq_min = 100 #Initial frequency
freq_max = 1000 #Final frequency
freq_step = 20 #Step between frequencies
freq = freq_min #Define the dynamical variable frequency

output_file_vpol = "vChanAntennaDef.dat" #Name of the output file
output_file_hpol = "hChanAntennaDef.dat"

# The name of the input file is broken into two pieces, since the actual path
# is dependent on the frequency, which changes as the files are looped over.
input_file_prefix_vpol = "antenna_vchan_"
input_file_prefix_hpol = "antenna_hchan_"
input_file_suffix = "_MHz.uan"
io_directory = "/Users/neutrino/Dropbox/GP_Antennas/OSU_internal/HornFreqSweep"


patt = open(output_file_vchan, "w")# Open file to write on it

while (freq <= freq_max): 
    input_file_vchan = (io_directory + input_file_prefix_vchan +
                       str(freq) + input_file_suffix) 

    #Prints the name of the input file
    print(input_file_vchan) 

    #Header for each iteration
    patt.write("freq : "+str(freq)+".00 MHz"+'\n'+"SWR : 1.965000"+'\n')
   
    if(freq==freq_min):
        
        # Write the header to the output file
        # Where we have a standard xyz-coordinate system, the
        # vertically polarized antenna will have its theta be
        # from the x-hat vector to the z-hat vector, and its
        # phi from the x-hat vector to the y-hat vector.
        # The horizontally polarized antenna will have it's
        # theta from the x-hat vector to the negative y-hat
        # vector, and its phi from the x-hat vector to the
        # z-hat vector.
        patt.write("Theta" + '\t' +
                   "Phi" + '\t' +
                   "Theta-pol Gain(dB)" + '\t' +
                   "Theta-pol Gain" + '\t' +
                   "Theta-pol Phase(deg)" + '\t' +
                   "Phi-pol Gain(dB)" + '\t' +
                   "Phi-pol Gain" + '\t' +
                   "Phi-pol Phase(deg)" + '\t' +
                   '\n')
   
    with open(input_file_vchan, 'r') as f:
        for _ in xrange(17): 
            next(f) #Skips the first 17 lines of the XF file,
                    #which are useless for us.
        for line in f:
                line = line.strip() #Identify lines.
                columns = line.split() #Identify columns.
                
                #Here I identify the veriables with their respective columns in
                #the XF file
                theta     = int(columns[0])
                phi       = int(columns[1]) 
                dB_t_gain = float(columns[2])
                t_gain    = math.pow(10, (dB_t_gain/10))
                dB_p_gain = float(columns[3])
                p_gain    = math.pow(10 , (dB_p_gain/10))
                t_phase   = float(columns[4])
                p_phase   = float(columns[5])
                #gain_tot = float(t_gain + p_gain)

                patt.write(str(theta) + '\t' +
                           str(phi) + '\t' +
                           str(dB_t_gain) + '\t' +
                           str(t_gain) + '\t' +
                           str(t_phase) + '\t' +
                           str(dB_p_gain) + '\t' +
                           str(p_gain) + '\t' +
                           str(p_phase) + '\n') #Write data on the new file
    freq += freq_step

#Close everything
f.close()
patt.close()
        
# Now redo that for all of the h-pol antenna files 

patt = open(output_file_hchan, "w")# Open file to write on it

while (freq <= freq_max): 
    input_file_hchan = (io_directory + input_file_prefix_hchan +
                       str(freq) + input_file_suffix) 

    #Prints the name of the input file
    print(input_file_hchan) 

    #Header for each iteration
    patt.write("freq : "+str(freq)+".00 MHz"+'\n'+"SWR : 1.965000"+'\n')
   
    if(freq==freq_min):
        
        # Write the header to the output file
        # Where we have a standard xyz-coordinate system, the
        # vertically polarized antenna will have it's theta be
        # from the x-hat vector to the z-hat vector, and it's
        # phi from the x-hat vector to the y-hat vector.
        # The horizontally polarized antenna will have it's
        # theta from the x-hat vector to the negative y-hat
        # vector, and it's phi from the x-hat vector to the
        # z-hat vector.
        patt.write("Theta" + '\t' +
                   "Phi" + '\t' +
                   "Theta-pol Gain(dB)" + '\t' +
                   "Theta-pol Gain" + '\t' +
                   "Theta-pol Phase(deg)" + '\t' +
                   "Phi-pol Gain(dB)" + '\t' +
                   "Phi-pol Gain" + '\t' +
                   "Phi-pol Phase(deg)" + '\t' +
                   '\n')
   
    with open(input_file_hchan, 'r') as f:
        for _ in xrange(17): 
            next(f) #Skips the first 17 lines of the XF file,
                    #which are useless for us.
        for line in f:
                line = line.strip() #Identify lines.
                columns = line.split() #Identify columns.
                
                #Here I identify the veriables with their respective columns in
                #the XF file
                theta     = int(columns[0])
                phi       = int(columns[1]) 
                dB_t_gain = float(columns[2])
                t_gain    = math.pow(10, (dB_t_gain/10))
                dB_p_gain = float(columns[3])
                p_gain    = math.pow(10 , (dB_p_gain/10))
                t_phase   = float(columns[4])
                p_phase   = float(columns[5])
                #gain_tot = float(t_gain + p_gain)

                patt.write(str(theta) + '\t' +
                           str(phi) + '\t' +
                           str(dB_t_gain) + '\t' +
                           str(t_gain) + '\t' +
                           str(t_phase) + '\t' +
                           str(dB_p_gain) + '\t' +
                           str(p_gain) + '\t' +
                           str(p_phase) + '\n') #Write data on the new file
    freq += freq_step

#Close everything
f.close()
patt.close()
 
