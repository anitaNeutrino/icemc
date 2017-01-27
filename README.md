# [icemc](https://github.com/anitaNeutrino/icemc) is a tool to simulate neutrino interactions in the ice #

## Prerequisites ##

To compile icemc you (should) only need a working [ROOT](https://root.cern.ch/downloading-root) installation.

Optional extras include the ANITA libraries [libRootFftwWrapper](https://github.com/nichol77/libRootFftwWrapper) and [eventReaderRoot](https://github.com/anitaNeutrino/eventReaderRoot) libraries.
The recommended way to install these is with the [anitaBuildTool](https://github.com/anitaNeutrino/anitaBuildTool) which downloads and compiles the core ANITA software libraries.
(See the ANITA-3 notes section for more information.)

## Getting icemc ##

``` bash
git clone https://github.com/anitaNeutrino/icemc
```
``` bash
make
```
The Makefile is currently in a sorry state and may require some love by uncommenting or commenting out various flags.
Makefiles are a dark art so send an email, or message us in the slack [#anita_simulation](https://anitamission.slack.com/messages/anita_simulation/) or [#github_commits](https://anitamission.slack.com/messages/github_commits/) rooms if you need help.


## Running icemc ##

From the command line do
``` bash
    ./icemc -i {inputFile} -o {outputDirectory} -r {runNumber} -n {numberOfNeutrinos} -t {triggerThreshold} -e {energyExponent}
```
If parameters are not specified

   * inputs from inputs.conf are used
   * the output directory is output
   * the run number is
   * number of neutrinos to generate as is input file
   * trigger threshold for full band as defined in anita.hh for that flight

## ANITA formatted data output ##

To produce Anita-like output files, you must have the [libRootFftwWrapper](https://github.com/nichol77/libRootFftwWrapper) and [eventReaderRoot](https://github.com/anitaNeutrino/eventReaderRoot) libraries installed before compiling icemc.

After that you will need to find and uncomment the line in icemc.cc `#define ANITA3_EVENTREADER`.
The next time you compile and run icemc, the ANITA-like outputs will be produced automatically.
The outputs also include the TruthAnitaEvent tree with the Monte Carlo truth information about the simulated events.

## ANITA-3 notes ##

 * If you want to generate anita-3 simulations with the impulse response and noise you should turn on the flags in the input file as by default they are off

* If you want to generate a pure thermal noise sample you should to turn on these flags:
     * Set signal to 0 to measure noise hits
     * skip cuts on neutrinos that won’t reach the payload (this is just to have it running quicker as we 0’ed the signal anyway)
     * Min bias flag (so that everything gets saved)
