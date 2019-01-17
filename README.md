# icemc
[icemc](https://github.com/anitaNeutrino/icemc) is a tool to simulate neutrino interactions in the ice.

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
Currently the Makefile tries to figure out whether to compile against the ANITA libraries by looking for the `ANITA_UTIL_INSTALL_DIR` environment variable.
The Makefile may require some love by uncommenting or commenting out various site specific flags at the top.
Makefiles are a dark art so if you need help send an email or message someone in the slack [#anita_simulation](https://anitamission.slack.com/messages/anita_simulation/) or [#github_commits](https://anitamission.slack.com/messages/github_commits/) rooms.

## Running icemc ##

To run icemc you need first to define two environment variables and then add them to your (DY)LD_LIBRARY_PATH, so that icemc knows where to look for its input files:
 
   * ICEMC_SRC_DIR should point to the directory where all your source code is (i.e. this directory)
   * ICEMC_BUILD_DIR should point to the directory where your exectutables and .pcm are (in case you run icemc by itself it's again this directory, but in case you installed icemc from the anitaBuildTool, it should be something different)

An example is:
```bash
export ICEMC_SRC_DIR=/path/to/anitaBuildTool/components/icemc/
export ICEMC_BUILD_DIR=/path/to/anitaBuildTool/build/components/icemc/
export DYLD_LIBRARY_PATH=${ICEMC_SRC_DIR}:${ICEMC_BUILD_DIR}:${DYLD_LIBRARY_PATH}
```

If you use the anitaBuildTool, it's a good practice to add these lines to your setup script.


To run icemc do:
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

After that you will need to find and uncomment the site specific flag `ANITA3_EVENTREADER=1` in the Makefile.
The next time you compile and run icemc, the ANITA-like outputs will be produced automatically.
The outputs also include the TruthAnitaEvent tree with the Monte Carlo truth information about the simulated events.

## ANITA-3 notes ##

 * If you want to generate anita-3 simulations with the impulse response and noise you should turn on the flags in the input file as by default they are off

* If you want to generate a pure thermal noise sample you should to turn on these flags:
     * Set signal to 0 to measure noise hits
     * skip cuts on neutrinos that won’t reach the payload (this is just to have it running quicker as we 0’ed the signal anyway)
     * Min bias flag (so that everything gets saved)
