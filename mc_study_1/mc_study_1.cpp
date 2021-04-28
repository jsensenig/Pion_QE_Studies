//
// Created by Jon Sensenig on 4/23/21.
//

#include "mc_study_1.h"
#include "../utilities/Utilities.hpp"
#include "TTree.h"
#include "TVector3.h"
#include <iostream>


void run_mc_study( std::string in_file, Histograms &hists ) {

  TFile *proc_file = TFile::Open( in_file.c_str() );

  if( !proc_file -> IsOpen() ) {
    std::cout << "File " << in_file << " not open!" << std::endl;
    return;
  }

  TTree* tree = (TTree*)proc_file -> Get("pionana/beamana");

  int true_beam_PDG;
  double true_beam_endX;
  double true_beam_endY;
  double true_beam_endZ;
  double true_beam_endP;
  double reco_beam_endX;
  double reco_beam_endY;
  double reco_beam_endZ;
  double reco_beam_startX;
  double reco_beam_startY;
  double reco_beam_startZ;
  double reco_beam_len;
  int reco_beam_type;
  int true_daughter_nNeutron;
  int true_daughter_nProton;
  int true_daughter_nPiMinus;
  int true_daughter_nPiPlus;
  int true_daughter_nPi0;
  double true_beam_endPx, true_beam_endPy, true_beam_endPz;
  std::string *true_beam_endProcess = new std::string;
  std::vector<double> *true_beam_daughter_startPx = new std::vector<double>;
  std::vector<double> *true_beam_daughter_startPy = new std::vector<double>;
  std::vector<double> *true_beam_daughter_startPz = new std::vector<double>;
  std::vector<int> *true_beam_daughter_PDG = new std::vector<int>;


  tree->SetBranchAddress("true_beam_PDG", &true_beam_PDG);
  tree->SetBranchAddress("true_beam_endX", &true_beam_endX);
  tree->SetBranchAddress("true_beam_endY", &true_beam_endY);
  tree->SetBranchAddress("true_beam_endZ", &true_beam_endZ);
  tree->SetBranchAddress("reco_beam_endX", &reco_beam_endX);
  tree->SetBranchAddress("reco_beam_endY", &reco_beam_endY);
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  tree->SetBranchAddress("reco_beam_len", &reco_beam_len);
  tree->SetBranchAddress("reco_beam_type", &reco_beam_type);
  tree->SetBranchAddress("reco_beam_startX", &reco_beam_startX);
  tree->SetBranchAddress("reco_beam_startY", &reco_beam_startY);
  tree->SetBranchAddress("reco_beam_startZ", &reco_beam_startZ);
  tree->SetBranchAddress("true_daughter_nNeutron", &true_daughter_nNeutron);
  tree->SetBranchAddress("true_daughter_nProton", &true_daughter_nProton);
  tree->SetBranchAddress("true_daughter_nPiMinus", &true_daughter_nPiMinus);
  tree->SetBranchAddress("true_daughter_nPiPlus", &true_daughter_nPiPlus);
  tree->SetBranchAddress("true_daughter_nPi0", &true_daughter_nPi0);
  tree->SetBranchAddress("true_beam_endProcess", &true_beam_endProcess);
  tree->SetBranchAddress("true_beam_endPx", &true_beam_endPx);
  tree->SetBranchAddress("true_beam_endPy", &true_beam_endPy);
  tree->SetBranchAddress("true_beam_endPz", &true_beam_endPz);
  tree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);
  tree->SetBranchAddress("true_beam_daughter_startPx", &true_beam_daughter_startPx);
  tree->SetBranchAddress("true_beam_daughter_startPy", &true_beam_daughter_startPy);
  tree->SetBranchAddress("true_beam_daughter_startPz", &true_beam_daughter_startPz);
  tree->SetBranchAddress("true_beam_endP", &true_beam_endP);


  size_t nevts = tree -> GetEntries();
  std::cout << "Processing " << nevts << " events." << std::endl;

  for ( size_t i = 0; i < nevts; i++ ) {
    tree->GetEntry( i );

    double pi_ke = -1;

    // Check what is included in the definition of the "pi+inelastic"
    if( true_beam_endProcess->compare("pi+Inelastic") == 0 ) {
      hists.th1_hists["hPrimaryPdg"] -> Fill( true_beam_PDG );
      hists.th1_hists["hDaughterPiP"] -> Fill( true_daughter_nPiPlus );
      hists.th1_hists["hDaughterPiM"] -> Fill( true_daughter_nPiMinus );
      hists.th1_hists["hDaughterPi0"] -> Fill( true_daughter_nPi0 );
      hists.th1_hists["hDaughterP"] -> Fill( true_daughter_nProton );
      hists.th1_hists["hDaughterN"] -> Fill( true_daughter_nNeutron );

      // All in-elastic pi events
      pi_ke = utils::CalculateKE( true_beam_endP*1.e3, utils::pdg::pdg2mass(utils::pdg::kPdgPiP) );
      hists.th1_hists["hPiInElKe"] -> Fill( pi_ke );

    }

    bool true_qe = true_beam_endProcess->compare("pi+Inelastic") == 0 &&
                   true_daughter_nPiPlus == 1 &&
                   true_daughter_nPiMinus == 0 &&
                   true_daughter_nPi0 == 0 &&
                   ( true_daughter_nNeutron > 0 || true_daughter_nProton > 0 );

    // Only look at true QE events
    if( !true_qe ) continue;

    // Nucleon multiplicity
    hists.th2_hists["hnProtonNeutron"] -> Fill( true_daughter_nProton, true_daughter_nNeutron );

    // Scattering angle of the daughter pion
    int pi_idx = utils::FindIndex( *true_beam_daughter_PDG, utils::pdg::kPdgPiP );
    double pi_angle = scatter_angle( true_beam_endPx, true_beam_endPy, true_beam_endPz,
                      true_beam_daughter_startPx->at(pi_idx), true_beam_daughter_startPy->at(pi_idx),
                      true_beam_daughter_startPz->at(pi_idx));
    hists.th1_hists["hPiAngle"] -> Fill( pi_angle );
    hists.th1_hists["hPiCosAngle"] -> Fill( TMath::Cos(pi_angle) );

    // Only pi QE events
    hists.th1_hists["hPiQeKe"] -> Fill( pi_ke );
    hists.th1_hists["hFracPiQeKe"] -> Fill( pi_ke );

  }

  // Divide pi QE / all pi In-elastic histograms to get fraction of pion QE events
  hists.th1_hists["hFracPiQeKe"] -> Divide( hists.th1_hists["hPiInElKe"].get() );

  // Clean up
  delete true_beam_endProcess;
  delete true_beam_daughter_PDG;
  delete true_beam_daughter_startPx;
  delete true_beam_daughter_startPy;
  delete true_beam_daughter_startPz;
  proc_file -> Close();

}

double scatter_angle( double x1, double y1, double z1, double x2, double y2, double z2 ) {

  TVector3 in( x1, y1, z1 );
  TVector3 out( x2, y2, z2 );
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

  std::cout << "Starting MC study!" << std::endl;
  run_mc_study( input_file, hists );

  std::cout << "Writing histograms to " << output_file << std::endl;
  // Write histograms ot file
  hists.WriteHistos( output_file );

  return 0;

}