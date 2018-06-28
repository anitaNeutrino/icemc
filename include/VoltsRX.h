#ifndef VOLTS_RX_H
#define VOLTS_RX_H


#include "TObject.h"
#include "anita.hh"

namespace icemc {

  /**
   * @class VoltsRX
   * @brief Voltage seen at the antennas
   */
  class VoltsRX {
  public:
    VoltsRX(int nRX);
    void reset();
    
    double max;		///< max voltage seen on an antenna - just for debugging purposes
    double ave;		///< ave voltage seen on an antenna,  among hit antennas
    double sum;		///< ave voltage seen on an antenna,  among hit antennas

    double max_highband;	///< max voltage seen on an antenna - just for debugging purposes
    double max_lowband;	///< max voltage seen on an antenna - just for debugging purposes

    const int fNRX;
    std::vector<std::array<double, Anita::HALFNFOUR> > rfcm_lab_e_all;
    std::vector<std::array<double, Anita::HALFNFOUR> > rfcm_lab_h_all;
    
    ClassDef(VoltsRX, 2);
  };
}

#endif
