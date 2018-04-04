#ifndef VOLTS_RX_H
#define VOLTS_RX_H


#include "TObject.h"

namespace icemc {

  /**
   * @class VoltsRX
   * @brief Voltage seen at the antennas
   */
  class VoltsRX {
  public:
    VoltsRX();
    
    double max;		///< max voltage seen on an antenna - just for debugging purposes
    double ave;		///< ave voltage seen on an antenna,  among hit antennas
    double sum;		///< ave voltage seen on an antenna,  among hit antennas

    double max_highband;	///< max voltage seen on an antenna - just for debugging purposes
    double max_lowband;	///< max voltage seen on an antenna - just for debugging purposes

    static const int nRX = 48;
    static const int nSamp = 512;
    
    double rfcm_lab_e_all[nRX][nSamp];
    double rfcm_lab_h_all[nRX][nSamp];

    ClassDef(VoltsRX, 1);
  };
}

#endif
