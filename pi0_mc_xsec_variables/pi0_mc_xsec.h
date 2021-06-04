//
// Created by Jon Sensenig on 4/23/21.
//

#ifndef PI0_MC_STUDY_H
#define PI0_MC_STUDY_H

#include "../utilities/Histograms.hpp"

/// Pi0 Data Structure
struct pi0_proxy {
  std::pair<int,int> id;            // Gamma ID
  std::pair<double,double> energy;  // Gamma E [MeV]
  std::pair<double,double> angle;   // Gamma angle [rad]
  double open_angle;                // Gamma open angle [rad]

  void reset() { // Reset all variables
    id         = {0, 0};
    energy     = {0., 0.};
    angle      = {0., 0.};
    open_angle = 0.;
  }
};

void run_pi0_mc_xsec( std::string in_file, Histograms &hists );

std::map<int, double> daughter_pi0_energy( Histograms& hists );

void extract_xsec( pi0_proxy& pi0 );

double pi0_angle( pi0_proxy& pi0 );

double pi0_mom( pi0_proxy& pi0 );

double pi0_energy( pi0_proxy& pi0 );

void plot_single_pi0( pi0_proxy& pi0, Histograms& hists, size_t idx );

void plot_all_pi0( pi0_proxy& pi0, double pi0_energy, Histograms& hists, size_t idx );

double open_angle( double px1, double py1, double pz1, double px2, double py2, double pz2 );

void extract_xsec( size_t nevts );

void plot_xsec( std::vector<std::vector<double>> &xsec, std::vector<double> &angle, std::vector<double> &energy );

void clean_pointers();

/// Histograms for X-section
auto energy_hist = std::make_unique<TH1D>( "pi0_energy", "Pi0 Energy;E_pi0 [MeV];Count", 12, 0., 1200. );
auto angle_hist = std::make_unique<TH1D>( "pi0_angle", "Pi0 Angle;#theta_pi0 [deg];Count", 18, 0., 180. );
auto energy_angle_hist = std::make_unique<TH2D>( "pi0_energy_angle", "Pi0 Angle vs Energy;E [MeV];#theta_pi0 [deg]", 20, 0., 2000., 9, 0., 90. );

/// Data Variables
int true_daughter_nPi0;
int true_daughter_nNeutron;
int true_daughter_nProton;
int true_daughter_nPiMinus;
int true_daughter_nPiPlus;
int true_beam_PDG;
double true_beam_endP;
double true_beam_endX, true_beam_endY, true_beam_endZ;
double true_beam_endPx, true_beam_endPy, true_beam_endPz;

auto true_beam_daughter_PDG = new std::vector<int>;
auto true_beam_daughter_ID = new std::vector<int>;
auto true_beam_daughter_startP = new std::vector<double>;
auto true_beam_daughter_startPx = new std::vector<double>;
auto true_beam_daughter_startPy = new std::vector<double>;
auto true_beam_daughter_startPz = new std::vector<double>;
auto true_beam_Pi0_decay_PDG = new std::vector<int>;
auto true_beam_Pi0_decay_parID = new std::vector<int>;
auto true_beam_Pi0_decay_startP = new std::vector<double>;
auto true_beam_Pi0_decay_startPx = new std::vector<double>;
auto true_beam_Pi0_decay_startPy = new std::vector<double>;
auto true_beam_Pi0_decay_startPz = new std::vector<double>;
auto true_beam_Pi0_decay_startX = new std::vector<double>;
auto true_beam_Pi0_decay_startY = new std::vector<double>;
auto true_beam_Pi0_decay_startZ = new std::vector<double>;
auto true_beam_Pi0_decay_len = new std::vector<double>;

auto true_beam_endProcess = new std::string;

#endif //PI0_MC_STUDY_H

