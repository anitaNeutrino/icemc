////////////
// rx.hpp //
////////////

#ifndef RX_HPP
#define RX_HPP 1

// Standard Library #includes
#include <vector>

// ROOT Library #includes
#include "TObject.h"
//#include "TCint.h"

class RX : public TObject {
public:
    unsigned phi_sector;
    unsigned layer;
    double x, y, z;
    std::vector <double>* waveform;
    std::vector <double>* digitized;
	RX (void) : waveform (NULL), digitized (NULL) {}
    ~RX (void) {delete waveform; delete digitized;}
    ClassDef(RX, 1);
};

#endif
