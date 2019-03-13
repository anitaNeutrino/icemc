#ifndef FAVA_ENTRY_H
#define  FAVA_ENTRY_H

#include "TObjString.h"

struct FAVAEntry {
    int week; 
    int num; 
    double ra; 
    double dec;
    double r95; 
    int met_tmin;
    int met_tmax;
    int unix_tmin;
    int unix_tmax;
    int le_nev;
    int he_nev;
    double he_avnev;
    double le_avnev;
    double he_ts;
    double le_ts;
    double he_ts_sigma;
    double le_ts_sigma;
    double he_sigma;
    double le_sigma;
    double le_flux;
    double he_flux;
    double le_flux_err;
    double he_flux_err;
    double le_index;
    double he_index;
    double z; 
    bool tevcat; 
    TObjString association;
    TObjString association_3fgl;
    TObjString source_class;
    };

#endif 


