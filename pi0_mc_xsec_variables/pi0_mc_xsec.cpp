//
// Created by Jon Sensenig on 4/23/21.
//

#include "pi0_mc_xsec.h"
#include "../utilities/Utilities.hpp"
#include "TTree.h"
#include "TVector3.h"
#include <iostream>


void run_pi0_mc_xsec( std::string in_file, Histograms &hists ) {

  TFile *proc_file = TFile::Open( in_file.c_str() );

  if( !proc_file -> IsOpen() ) {
    std::cout << "File " << in_file << " not open!" << std::endl;
    return;
  }

  TTree* tree = (TTree*)proc_file -> Get("pionana/beamana");

  int true_daughter_nPi0;
  int true_daughter_nNeutron;
  int true_daughter_nProton;
  int true_daughter_nPiMinus;
  int true_daughter_nPiPlus;
  int true_beam_PDG;
  double true_beam_endP;
  double true_beam_endX, true_beam_endY, true_beam_endZ;
  double true_beam_endPx, true_beam_endPy, true_beam_endPz;

  std::vector<int> *true_beam_daughter_PDG = new std::vector<int>;
  std::vector<int> *true_beam_daughter_ID = new std::vector<int>;
  std::vector<double> *true_beam_daughter_startP = new std::vector<double>;
  std::vector<double> *true_beam_daughter_startPx = new std::vector<double>;
  std::vector<double> *true_beam_daughter_startPy = new std::vector<double>;
  std::vector<double> *true_beam_daughter_startPz = new std::vector<double>;
  std::vector<int> *true_beam_Pi0_decay_PDG = new std::vector<int>;
  std::vector<int> *true_beam_Pi0_decay_parID = new std::vector<int>;
  std::vector<double> *true_beam_Pi0_decay_startP = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_startPx = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_startPy = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_startPz = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_startX = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_startY = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_startZ = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_len = new std::vector<double>;

  std::string *true_beam_endProcess = new std::string;

  tree->SetBranchAddress("true_daughter_nPi0", &true_daughter_nPi0);
  tree->SetBranchAddress("true_daughter_nNeutron", &true_daughter_nNeutron);
  tree->SetBranchAddress("true_daughter_nProton", &true_daughter_nProton);
  tree->SetBranchAddress("true_daughter_nPiMinus", &true_daughter_nPiMinus);
  tree->SetBranchAddress("true_daughter_nPiPlus", &true_daughter_nPiPlus);
  tree->SetBranchAddress("true_beam_PDG", &true_beam_PDG);
  tree->SetBranchAddress("true_beam_endP", &true_beam_endP);
  tree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);
  tree->SetBranchAddress("true_beam_daughter_ID", &true_beam_daughter_ID);
  tree->SetBranchAddress("true_beam_daughter_startP", &true_beam_daughter_startP);
  tree->SetBranchAddress("true_beam_endProcess", &true_beam_endProcess);
  tree->SetBranchAddress("true_beam_endX", &true_beam_endX);
  tree->SetBranchAddress("true_beam_endY", &true_beam_endY);
  tree->SetBranchAddress("true_beam_endZ", &true_beam_endZ);
  tree->SetBranchAddress("true_beam_endPx", &true_beam_endPx);
  tree->SetBranchAddress("true_beam_endPy", &true_beam_endPy);
  tree->SetBranchAddress("true_beam_endPz", &true_beam_endPz);
  tree->SetBranchAddress("true_beam_daughter_startPx", &true_beam_daughter_startPx);
  tree->SetBranchAddress("true_beam_daughter_startPy", &true_beam_daughter_startPy);
  tree->SetBranchAddress("true_beam_daughter_startPz", &true_beam_daughter_startPz);

  tree->SetBranchAddress("true_beam_Pi0_decay_PDG", &true_beam_Pi0_decay_PDG);
  tree->SetBranchAddress("true_beam_Pi0_decay_startP", &true_beam_Pi0_decay_startP);
  tree->SetBranchAddress("true_beam_Pi0_decay_parID", &true_beam_Pi0_decay_parID);
  tree->SetBranchAddress("true_beam_Pi0_decay_startPx", &true_beam_Pi0_decay_startPx);
  tree->SetBranchAddress("true_beam_Pi0_decay_startPy", &true_beam_Pi0_decay_startPy);
  tree->SetBranchAddress("true_beam_Pi0_decay_startPz", &true_beam_Pi0_decay_startPz);
  tree->SetBranchAddress("true_beam_Pi0_decay_startX", &true_beam_Pi0_decay_startX);
  tree->SetBranchAddress("true_beam_Pi0_decay_startY", &true_beam_Pi0_decay_startY);
  tree->SetBranchAddress("true_beam_Pi0_decay_startZ", &true_beam_Pi0_decay_startZ);
  tree->SetBranchAddress("true_beam_Pi0_decay_len", &true_beam_Pi0_decay_len);

  size_t nevts = tree -> GetEntries();
  std::cout << "Processing " << nevts << " events." << std::endl;

  for ( size_t i = 0; i < nevts; i++ ) {
    tree->GetEntry( i );

    double pi0_mass = utils::pdg::pdg2mass( utils::pdg::kPdgPi0 );
    double pi_mass = utils::pdg::pdg2mass( utils::pdg::kPdgPiP );
    int pi0_idx = utils::FindIndex<int>( *true_beam_daughter_PDG, utils::pdg::kPdgPi0 );

    // Define true CEX
    bool true_cex = true_beam_endProcess->compare("pi+Inelastic") == 0 &&
                   true_daughter_nPi0 > 0 && true_daughter_nPiMinus == 0 &&
                   true_daughter_nPiPlus == 0;

    if( !true_cex ) continue;

    std::map<int,double> pi0_energy_map;
    for( size_t i = 0; i < true_beam_daughter_PDG->size(); i++ ) {

      // Only look at pi0
      if( true_beam_daughter_PDG->at(i) != utils::pdg::kPdgPi0 ) continue;
      // Plot the daughter pi0
      double daughter_pi0_energy = utils::CalculateE( true_beam_daughter_startP->at( i ) * 1.e3, pi0_mass );
      hists.th2_hists["hPiPPi0P"] -> Fill( true_beam_endP*1.e3, true_beam_daughter_startP->at( i )*1.e3 );
      if( true_daughter_nPi0 == 1 ) {
        hists.th2_hists["hPiP1Pi0P"] -> Fill( true_beam_endP*1.e3, true_beam_daughter_startP->at( i )*1.e3 );
      }
      hists.th1_hists["hPi0E"]->Fill( daughter_pi0_energy );
      pi0_energy_map[true_beam_daughter_ID->at(i)] = daughter_pi0_energy;
    }

    for( size_t i = 0; i < true_beam_Pi0_decay_parID->size()/2; i++ ) {
      size_t idx = i * 2; // Gamma from the same pi0 are adjacent in index, i.e., i and i+1
      if( true_beam_Pi0_decay_parID -> at( idx ) != true_beam_Pi0_decay_parID -> at( idx+1 ) ) {
        std::cout << "Something's wrong, decay gamma IDs do not match!" << std::endl;
        std::cout << "PDG 1/2 " << true_beam_Pi0_decay_PDG->at(idx) << "/" << true_beam_Pi0_decay_PDG->at(idx+1) << std::endl;
      }
      // Decay gammas momentum and opening angle
      double gamma1 = true_beam_Pi0_decay_startP -> at( idx ) * 1.e3;
      double gamma2 = true_beam_Pi0_decay_startP -> at( idx+1 ) * 1.e3;
      double angle = open_angle( true_beam_Pi0_decay_startPx->at(idx), true_beam_Pi0_decay_startPy->at(idx),
                                 true_beam_Pi0_decay_startPz->at(idx), true_beam_Pi0_decay_startPx->at(idx+1),
                                 true_beam_Pi0_decay_startPy->at(idx+1), true_beam_Pi0_decay_startPz->at(idx+1) );

      double beam_end_pos = utils::Distance(true_beam_endX, true_beam_endY, true_beam_endZ);
      double gamma_pos1 = utils::Distance(true_beam_Pi0_decay_startX->at(idx), true_beam_Pi0_decay_startY->at(idx), true_beam_Pi0_decay_startZ->at(idx));
      double gamma_pos2 = utils::Distance(true_beam_Pi0_decay_startX->at(idx+1), true_beam_Pi0_decay_startY->at(idx+1), true_beam_Pi0_decay_startZ->at(idx+1));

      hists.th1_hists["hGammaOpenAngle"] -> Fill( angle );
      hists.th2_hists["hPi0EGammaOpenAngle"] -> Fill( TMath::RadToDeg()*angle, pi0_energy_map.at(true_beam_Pi0_decay_parID->at( idx )) );
      hists.th1_hists["hLeadGammaP"] -> Fill( std::max( gamma1, gamma2 ) );
      hists.th1_hists["hSubLeadGammaP"] -> Fill( std::min( gamma1, gamma2 ) );
      hists.th2_hists["hGammaP"] -> Fill( std::max( gamma1, gamma2 ), std::min( gamma1, gamma2 ) );
      hists.th2_hists["hGammaR"] -> Fill( gamma_pos1-beam_end_pos, gamma_pos2-beam_end_pos );
      hists.th2_hists["hGammaLen"] -> Fill( true_beam_Pi0_decay_len->at(idx), true_beam_Pi0_decay_len->at(idx+1) );

      if( true_daughter_nPi0 > 1 ) continue;

      // Angular momentum calculation
      double angle_deg = angle * TMath::RadToDeg();
      // Momentum from polynomial fit Ref https://arxiv.org/pdf/1511.00941.pdf
      double p_poly = 2202.3 - 94.9*angle_deg + 2.1*pow(angle_deg, 2) - 0.025*pow(angle_deg, 3) + 0.00017*pow(angle_deg, 4)
               - 6.0e-7*pow(angle_deg, 5) + 8.5e-10*pow(angle_deg, 6);
      double p_root = sqrt( pow(gamma1,2) + pow(gamma2,2) + 2*gamma1*gamma2*cos(angle) );

      // Energy true and calculated
      double pi0_true_energy = utils::CalculateE( true_beam_daughter_startP->at( pi0_idx ) * 1.e3, pi0_mass );
      double pi0_calc_energy = sqrt( pow(pi0_mass, 2) + pow(p_root, 2) );

      hists.th1_hists["hPolyPi0PError"] -> Fill( p_poly / (true_beam_daughter_startP->at( pi0_idx )*1.e3) );
      hists.th1_hists["hRootPi0PError"] -> Fill( p_root / (true_beam_daughter_startP->at( pi0_idx )*1.e3) );
      hists.th2_hists["hPi0TrueCalcEnergy"] -> Fill( pi0_calc_energy, pi0_true_energy );

      // Gamma angle wrt pi0
      double gamma1_open_angle = open_angle( true_beam_daughter_startPx->at(pi0_idx), true_beam_daughter_startPy->at(pi0_idx),
                                        true_beam_daughter_startPz->at(pi0_idx), true_beam_Pi0_decay_startPx->at(idx),
                                        true_beam_Pi0_decay_startPy->at(idx), true_beam_Pi0_decay_startPz->at(idx) );
      double gamma2_open_angle = open_angle( true_beam_daughter_startPx->at(pi0_idx), true_beam_daughter_startPy->at(pi0_idx),
                                        true_beam_daughter_startPz->at(pi0_idx), true_beam_Pi0_decay_startPx->at(idx+1),
                                        true_beam_Pi0_decay_startPy->at(idx+1), true_beam_Pi0_decay_startPz->at(idx+1) );

      // Gamma angle wrt incoming pi+
      double gamma1_angle = open_angle( true_beam_endPx, true_beam_endPy, true_beam_endPz,
                                             true_beam_Pi0_decay_startPx->at(idx), true_beam_Pi0_decay_startPy->at(idx),
                                             true_beam_Pi0_decay_startPz->at(idx) );
      double gamma2_angle = open_angle( true_beam_endPx, true_beam_endPy, true_beam_endPz,
                                             true_beam_Pi0_decay_startPx->at(idx+1), true_beam_Pi0_decay_startPy->at(idx+1),
                                             true_beam_Pi0_decay_startPz->at(idx+1) );

      // pi0 angle wrt incoming pi+
      double pi0_angle = open_angle( true_beam_endPx, true_beam_endPy, true_beam_endPz,
                                     true_beam_daughter_startPx->at(pi0_idx), true_beam_daughter_startPy->at(pi0_idx),
                                     true_beam_daughter_startPz->at(pi0_idx) );

      hists.th1_hists["hGammaDiff"] -> Fill( gamma1_open_angle - gamma2_open_angle );
      hists.th2_hists["hGammaPi0Angle"] -> Fill( gamma1_open_angle*TMath::RadToDeg(), gamma2_open_angle*TMath::RadToDeg() );

      // angles wrt incoming pi+
      double cos_pi0_angle_calc = ( gamma1 * cos(gamma1_angle) + gamma2 * cos(gamma2_angle) ) / p_root;
      double cos_pi0_angle_mc = cos( pi0_angle );
      hists.th2_hists["hPi0TrueCalc"] -> Fill( cos_pi0_angle_calc, cos_pi0_angle_mc );
      hists.th1_hists["hPi0CalcAngleDiff"] -> Fill( cos_pi0_angle_calc - cos_pi0_angle_mc );
    }

  }

  // Clean up
  delete true_beam_endProcess;
  delete true_beam_daughter_startP;
  delete true_beam_daughter_PDG;
  delete true_beam_daughter_ID;
  delete true_beam_daughter_startPx;
  delete true_beam_daughter_startPy;
  delete true_beam_daughter_startPz;
  delete true_beam_Pi0_decay_PDG;
  delete true_beam_Pi0_decay_startP;
  delete true_beam_Pi0_decay_parID;
  delete true_beam_Pi0_decay_startPx;
  delete true_beam_Pi0_decay_startPy;
  delete true_beam_Pi0_decay_startPz;
  delete true_beam_Pi0_decay_startX;
  delete true_beam_Pi0_decay_startY;
  delete true_beam_Pi0_decay_startZ;
  delete true_beam_Pi0_decay_len;
  proc_file -> Close();

}

double open_angle( double px1, double py1, double pz1, double px2, double py2, double pz2 ) {

  TVector3 in( px1, py1, pz1 );
  TVector3 out( px2, py2, pz2 );
  // Angle from dot product definition
  return in.Angle( out );
}

int main() {

  std::string input_file = "../../../pionana_Prod4_mc_1GeV_1_14_21.root";
  TString output_file = "out.root";
  std::string hists_config = "../hists.json";

  // Configure histograms
  Histograms hists;
  hists.ConfigureHistos( hists_config );

  std::cout << "Starting pi0 xsec study!" << std::endl;
  run_pi0_mc_xsec( input_file, hists );

  std::cout << "Writing histograms to " << output_file << std::endl;
  // Write histograms ot file
  hists.WriteHistos( output_file );

  return 0;

}
