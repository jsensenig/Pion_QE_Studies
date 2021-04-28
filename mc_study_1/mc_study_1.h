//
// Created by Jon Sensenig on 4/23/21.
//

#ifndef MC_STUDY_1_H
#define MC_STUDY_1_H

#include "../utilities/Histograms.hpp"
#include "TLorentzVector.h"

void run_mc_study( std::string in_file );
void calculate_energy_loss( TLorentzVector &k_in, TLorentzVector &k_out, Histograms &hists );
double scatter_angle( double x1, double y1, double z1, double x2, double y2, double z2 );


#endif //MC_STUDY_1_H

