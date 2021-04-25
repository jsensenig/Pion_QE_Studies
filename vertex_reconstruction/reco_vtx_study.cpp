//
// Created by Jon Sensenig on 4/23/21.
//

#include "reco_vtx_study.h"
#include "TTree.h"
#include <iostream>


void run_reco_vertex( std::string in_file, Histograms &hists ) {

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
  double reco_beam_endX;
  double reco_beam_endY;
  double reco_beam_endZ;

  tree->SetBranchAddress("true_beam_PDG", &true_beam_PDG);
  tree->SetBranchAddress("true_beam_endX", &true_beam_endX);
  tree->SetBranchAddress("true_beam_endY", &true_beam_endY);
  tree->SetBranchAddress("true_beam_endZ", &true_beam_endZ);
  tree->SetBranchAddress("reco_beam_endX", &reco_beam_endX);
  tree->SetBranchAddress("reco_beam_endY", &reco_beam_endY);
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);

  size_t nevts = tree -> GetEntries();
  std::cout << "Processing " << nevts << " events." << std::endl;

  for ( size_t i = 0; i < nevts; i++ ) {
    tree->GetEntry( i );

    double delta_x = true_beam_endX - reco_beam_endX;
    double delta_y = true_beam_endY - reco_beam_endY;
    double delta_z = true_beam_endZ - reco_beam_endZ;

    hists.th1_hists["hTrueBeamPdg"] -> Fill( true_beam_PDG );

    // Fill the histograms per particle code 
    switch ( abs(true_beam_PDG) ) {
      case utils::pdg::kPdgElectron :
        hists.th1_hists["hEndXDiffEl"] -> Fill( delta_x );
        hists.th1_hists["hEndYDiffEl"] -> Fill( delta_y );
        hists.th1_hists["hEndZDiffEl"] -> Fill( delta_z );
        break;
      case utils::pdg::kPdgGamma :
        hists.th1_hists["hEndXDiffGa"] -> Fill( delta_x );
        hists.th1_hists["hEndYDiffGa"] -> Fill( delta_y );
        hists.th1_hists["hEndZDiffGa"] -> Fill( delta_z );
        break;
      case utils::pdg::kPdgMuon :
        hists.th1_hists["hEndXDiffMu"] -> Fill( delta_x );
        hists.th1_hists["hEndYDiffMu"] -> Fill( delta_y );
        hists.th1_hists["hEndZDiffMu"] -> Fill( delta_z );
        break;
      case utils::pdg::kPdgPiP :
        hists.th1_hists["hEndXDiffPi"] -> Fill( delta_x );
        hists.th1_hists["hEndYDiffPi"] -> Fill( delta_y );
        hists.th1_hists["hEndZDiffPi"] -> Fill( delta_z );
        break;
      case utils::pdg::kPdgKP :
        hists.th1_hists["hEndXDiffKa"] -> Fill( delta_x );
        hists.th1_hists["hEndYDiffKa"] -> Fill( delta_y );
        hists.th1_hists["hEndZDiffKa"] -> Fill( delta_z );
        break;
      case utils::pdg::kPdgProton :
        hists.th1_hists["hEndXDiffPro"] -> Fill( delta_x );
        hists.th1_hists["hEndYDiffPro"] -> Fill( delta_y );
        hists.th1_hists["hEndZDiffPro"] -> Fill( delta_z );
        break;
      default :
        std::cout << "Missed PDG " << true_beam_PDG << std::endl;
    }

  }

}

int main() {

  std::string input_file = "../../pionana_Prod4_mc_1GeV_1_14_21.root";
  TString output_file = "out.root";
  std::string hists_config = "hists.json";

  // Configure histograms
  Histograms hists;
  hists.ConfigureHistos( hists_config );

  std::cout << "Starting vertex reconstruction study!" << std::endl;
  run_reco_vertex( input_file, hists );

  std::cout << "Writing histograms to " << output_file << std::endl;
  // Write histograms ot file
  hists.WriteHistos( output_file );

  return 0;

}