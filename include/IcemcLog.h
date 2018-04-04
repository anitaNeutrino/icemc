#ifndef TEXT_OUTPUT_H
#define TEXT_OUTPUT_H

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

namespace icemc {

  /**
   * @class Log
   * @brief A pretty (and probably over-engineered) logger for icemc, 
   * which can potentially hold a number of text files.
   * 
   * This logger has the operator<< overloaded and should behave like an std::ostream, 
   * except that it will pass text to both a log file and std::cout or std::cerr as well, 
   * depending on whether the message is a warning/error.
   * 
   * One can log into the generic output and terminal with 
   * 
   *    log << "This text goes into both" << std::endl;
   * 
   * One can log directly into a file as usual (i.e. not appear in the terminal) with 
   * 
   *    log.foutput << "This just goes in the file" << std::endl;
   * 
   * One can pass flags for info, warning or error, which will change append the word
   * "Info!", "Warning!" or "Error!" as well as changing (only) the terminal color
   * In order to reset the color and stream (back to cout), pass std::endl at the end of the line.
   * 
   *    log << icemc::Log::error << "This text is for an error, will go to cerr and will be red" << std::endl;
   * 
   */

  class Log {
  public:

    enum severity {
      info,
      warning,
      error
    };

    /** 
     * Constructor
     * 
     * @param outputDir The directory in which all the contained logs will be written
     * @param run The current run, used for naming files
     */
    Log(const char* outputDir, int run);

    /** 
     * Destructor, makes sure the terminal colors are reset
     */
    virtual ~Log();

    /** 
     * Prints the output stream to std::cout and the default output file (foutput)
     * 
     * @param s is a streamable type
     */
    template <typename Streamable>
    Log& operator<<(const Streamable& s){      
      return message(s);
    }

    /** 
     * @brief Implements icemc::Log << std::endl, which resets the terminal colors
     * and sets the stream back to std::cout
     * 
     * How to make this behave more or less like a stream object with std::endl;
     * https://stackoverflow.com/questions/1134388/
     */
    typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
    typedef CoutType& (*StandardEndLine)(CoutType&);
    Log& operator<<(StandardEndLine manip){
      manip(getStream());

      if(fMustReset){	
	getStream() << getColorReset();
	fUseStdErr = false;
	fMustReset = false;
      }
      foutput << std::endl;
      return *this;
    }

    /** 
     * Turn on (or off) the text strings passed to cout/cerr which give color in ANSI terminals
     * 
     * @param useColorCodes true to have on, false to switch off, default is on (true)
     */
    void setUseColorCodes(bool useColorCodes = true){
      fUseColorCodes = useColorCodes;
    }

    /** 
     * Are the color codes turned on?
     * 
     * @return true if they are, false otherwise
     */
    bool getUseColorCodes() const { return fUseColorCodes;}    
    

    std::ofstream foutput;
    std::ofstream nu_out;
    std::ofstream veff_out;
    std::ofstream distanceout;
    std::ofstream outfile;
    std::ofstream forbrian;
    std::ofstream al_voltages_direct;
    std::ofstream eventsthatpassfile;
    std::ofstream fnumbers;
    std::ofstream fslac_viewangles;
    std::ofstream fslac_hitangles;

  private:

    /** 
     * Turns the corresponding cout stream a different color
     * 
     * @param s is info, warning, error
     * 
     * @return reference to self
     */
    Log& message(severity s);


    /** 
     * Pipes the message to the default output file and the default terminal stream
     * 
     * @param s a streamable template
     * 
     * @return reference to self
     */
    template<typename Streamable>
    Log& message(const Streamable& s) {      
      getStream() << s;
      foutput << s;
      return *this;
    }

    /** 
     * Get the terminal color reset code
     * @return the ANSI text string, if using color codes.
     */
    const char* getColorReset(){
      return fUseColorCodes ? "\x1b[0m" : "";
    }

    /** 
     * @brief Gets the correct terminal output stream.
     * 
     * Are we printing to std::cout or std::cerr?
     * @return cerr if fUseStdErr, cout otherwise
     */
    std::ostream& getStream(){
      if(!fUseStdErr){
	return std::cout;
      }
      else{
	return std::cerr;
      }
    }
    

    std::string fOutputDir;
    int fRun;
    bool fMustReset;
    bool fUseStdErr;
    bool fUseColorCodes;
  };
}


#endif
