#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TBranch.h>
#include <TH2.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TRandom.h>
#include <TEventList.h>
#include <TMultiLayerPerceptron.h>
#include <TComplex.h>
#include <TVirtualGeoPainter.h>
#include <TFile.h>
#include <TSystem.h>
#include <TClassTree.h>
#include <TPaveLabel.h>
#include <TCanvas.h>
#include<TGClient.h>
#include <RQ_OBJECT.h>
#include <TApplication.h>
#include <TRint.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <TF2.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <Riostream.h>
#include <TObjString.h>
#include <TChain.h>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TLatex.h"

void invMassMuMu()
  {

    TString infile="eicReconOutput/EICreconOut_JPsiMuMu_1ifb_10x250ep_Pruned.root";
    std::string filename = infile.Data();
    std::string beam_config;
    std::string marker = "_Pruned.root";

    size_t pos = filename.rfind(marker);
    beam_config = filename.substr(pos - 8, 8);

    std::cout << "Extracted beam config: " << beam_config << std::endl;

    std::string outfilename = "outputs/invMassOutput_" + beam_config + ".root";

    // Set output file for the histograms
    TFile *ofile = TFile::Open(outfilename.c_str(),"RECREATE");

    // Set up input file chain
    TChain *mychain = new TChain("events");
    mychain->Add(infile);

    // Initialize reader
    TTreeReader tree_reader(mychain);

    // Get Particle Information
    TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<double> partMomX(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> partMomY(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> partMomZ(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");

    // Get Reconstructed Track Information
    TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");

    // Get Associations Between MCParticles and ReconstructedChargedParticles
    TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "ReconstructedChargedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations.simID");

    // Define Eta Histograms

    TH1D *muonEta = new TH1D("muonEta", " #eta of Thrown Muons",100,-5.,5.);
    TH1D *matchedMuonEta = new TH1D("matchedMuonEta", " #eta of Thrown Muons That Have Matching Track",100,-5.,5.);
    TH1D *muonEff = new TH1D("muonEff","Efficency;  #eta",100,-5.,5.);
    
    TH1D *JPsiEta = new TH1D("JPsiEta", " #eta of Thrown J/PSi",100,-5.,5.);
    TH1D *matchedJPsiEta = new TH1D("matchedJPsiEta", " #eta of Thrown J/PSi That Have Matching Track",100,-5.,5.);
    TH1D *JPsiEff = new TH1D("JPsiEff","Efficency;  #eta",100,-5.,5.);

    // Define Momentum Histograms

    TH1D *muonMomHist = new TH1D("muonMomHist","Momentum of Thrown Muons;Momentum",100,-300.,300.);
    TH1D *matchedMuonMomHist = new TH1D("matchedMuonMomHist","Momentum of Thrown Muons That Have Matching Track",100,-300.,300.);
    TH1D *muonMomEff = new TH1D("muonMomEff","Efficency;Momentum",100,-300.,300.);

    TH1D *JPsiMomHist = new TH1D("JPsiMomHist","Momentum of Thrown J/PSi;Momentum",100,-300.,300.);
    TH1D *matchedJPsiMomHist = new TH1D("matchedJPsiMomHist","Momentum of Thrown J/PSi That Have Matching Track",100,-300.,300.);
    TH1D *JPsiMomEff = new TH1D("JPsiMomEff","Efficency;Momentum",100,-300.,300.);

    // Define 2D Histograms p vs eta

    TH2D *muonEtaMom = new TH2D("muonEtaMom","Thrown Muons;  #eta; Momentum",100,-5.,5.,100,-300.,300.);
    TH2D *matchedMuonEtaMom = new TH2D("matchedMuonEtaMom","Thrown Muons That Have Matching Track;  #eta; Momentum",100,-5.,5.,100,-300.,300.);

    TH2D *JPsiEtaMom = new TH2D("JPsiEtaMom","Thrown J/PSi;  #eta; Momentum",100,-5.,5.,100,-300.,300.);
    TH2D *matchedJPsiEtaMom = new TH2D("matchedJPsiEtaMom","Thrown J/PSi That Have Matching Track;  #eta; Momentum",100,-5.,5.,100,-300.,300.);

    // Define Delta R Histograms
    TH1D *matchedMuonTrackDeltaR = new TH1D("matchedMuonTrackDeltaR","Delta R Between Matching Thrown and Reconstructed Muon",5000,0.,5.);
    TH1D *matchedJPsiTrackDeltaR = new TH1D("matchedJPsiTrackDeltaR","Delta R Between Matching Thrown and Reconstructed J/PSi",5000,0.,5.);


    // Save 4Vectors of Muons
    int muonCount = 0;
    TLorentzVector trueMuPlus4Mom;
    TLorentzVector trueMuMinus4Mom;
    int reconstructedMuonCount = 0;
    TLorentzVector recMuPlus4Mom;
    TLorentzVector recMuMinus4Mom;

    TH1D *invMassMuMu = new TH1D("invMassMuMu","Invariant Mass of Muon Pairs",100,0.,5.); // Invariant mass histogram
    TH1D *matchedMuMuInvMass = new TH1D("matchedMuMuInvMass","Invariant Mass of Muon Pairs That Have Matching Track",100,0.,5.); // Invariant mass histogram

    while(tree_reader.Next()) { // Loop over events

      muonCount = 0; // Reset muon count
      trueMuPlus4Mom.SetPxPyPzE(0.,0.,0.,0.); // Reset 4-momenta
      trueMuMinus4Mom.SetPxPyPzE(0.,0.,0.,0.); // Reset 4-momenta
      reconstructedMuonCount = 0; // Reset reconstructed muon count
      recMuPlus4Mom.SetPxPyPzE(0.,0.,0.,0.); // Reset 4-momenta
      recMuMinus4Mom.SetPxPyPzE(0.,0.,0.,0.); // Reset 4-momenta 

      for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop over thrown particles
        {
          if(partGenStat[i] == 1) // Select stable thrown particles
            {
              int pdg = TMath::Abs(partPdg[i]);
  
              if(pdg == 13 || pdg == -13) // Look at muons
                {
                  TVector3 trueMom(partMomX[i],partMomY[i],partMomZ[i]);
                  TLorentzVector true4Mom(trueMom, sqrt(trueMom.Mag2() + 0.105658*0.105658)); // Muon mass = 0.105658 GeV/c^2

                  if (pdg == 13) // Positive muon
                    {
                      muonCount++; // Increment muon count
                      trueMuPlus4Mom += true4Mom; // Add to positive muon 4-momentum
                    }
                  else if (pdg == -13) // Negative muon
                    {
                      muonCount++; // Increment muon count
                      trueMuMinus4Mom += true4Mom; // Add to negative muon 4-momentum
                    }

                  float trueEta = trueMom.PseudoRapidity();
                  float truePhi = trueMom.Phi();
        
                  muonEta->Fill(trueEta);
                  muonMomHist->Fill(trueMom.Mag());
                  muonEtaMom->Fill(trueEta,trueMom.Mag());
                  // Loop over associations to find matching ReconstructedChargedParticle
                  for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
                    {
                      if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
                        {
                          TVector3 recMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle
                          TLorentzVector rec4Mom(recMom, sqrt(trueMom.Mag2() + 0.105658*0.105658)); // Muon mass = 0.105658 GeV/c^2

                          float recEta = recMom.PseudoRapidity();
                          float recPhi = recMom.Phi();

                          if (pdg == 13) // Positive muon
                            {
                              reconstructedMuonCount++; // Increment reconstructed muon count
                              recMuPlus4Mom += rec4Mom; // Add to positive muon 4-momentum
                            }
                          else if (pdg == -13) // Negative muon
                            {
                              reconstructedMuonCount++; // Increment reconstructed muon count
                              recMuMinus4Mom += rec4Mom; // Add to negative muon 4-momentum
                            }

                          // Check the distance between the thrown and reconstructed particle
                          float deltaEta = trueEta - recEta;
                          float deltaPhi = TVector2::Phi_mpi_pi(truePhi - recPhi);
                          float deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

                          matchedMuonTrackDeltaR->Fill(deltaR);

                          matchedMuonEta->Fill(recEta); // Plot the thrown eta if a matched ReconstructedChargedParticle was found
                          muonEff->Fill(recEta);
                          matchedMuonMomHist->Fill(recMom.Mag());
                          muonMomEff->Fill(recMom.Mag());

                          matchedMuonEtaMom->Fill(recEta,recMom.Mag());

                          matchedJPsiTrackDeltaR->Fill(deltaR); // Fill delta R histogram for J/PSi
                        }
                    }
                }       
            }
        }

      // Calculate invariant mass of muon pairs
      if (muonCount == 2) 
        {
          double invMass = (trueMuPlus4Mom + trueMuMinus4Mom).M(); // Invariant mass
          invMassMuMu->Fill(invMass); // Fill histogram

          JPsiEta->Fill(trueMuPlus4Mom.Eta());
          JPsiMomHist->Fill(trueMuPlus4Mom.P());
          JPsiEtaMom->Fill(trueMuPlus4Mom.Eta(), trueMuPlus4Mom.P());

          if (reconstructedMuonCount == 2)
            {
              double matchedInvMass = (recMuPlus4Mom + recMuMinus4Mom).M(); // Invariant mass
              matchedMuMuInvMass->Fill(matchedInvMass); // Fill histogram

              matchedJPsiEta->Fill(recMuPlus4Mom.Eta());
              JPsiEff->Fill(recMuPlus4Mom.Eta());
              matchedJPsiMomHist->Fill(recMuPlus4Mom.P());
              JPsiMomEff->Fill(recMuPlus4Mom.P());
              matchedJPsiEtaMom->Fill(recMuPlus4Mom.Eta(), recMuPlus4Mom.P());

            }
        }

    }

    // Calculate efficiency
    muonEff->Divide(muonEta);
    muonMomEff->Divide(muonMomHist);

    JPsiEff->Divide(JPsiEta);
    JPsiMomEff->Divide(JPsiMomHist);

    // Write histograms to output file
    ofile->cd();
    ofile->mkdir("muons");
    ofile->cd("muons");
    muonEta->Write();
    matchedMuonEta->Write();
    muonEff->Write();
    muonMomHist->Write();
    matchedMuonMomHist->Write();
    muonMomEff->Write();
    muonEtaMom->Write();
    matchedMuonEtaMom->Write();
    matchedMuonTrackDeltaR->Write();
    ofile->cd("..");
    ofile->mkdir("invariantMass");
    ofile->cd("invariantMass");
    invMassMuMu->Write();
    matchedMuMuInvMass->Write();
    ofile->cd("..");
    ofile->mkdir("JPsi");
    ofile->cd("JPsi");
    JPsiEta->Write();
    matchedJPsiEta->Write();
    JPsiEff->Write();
    JPsiMomHist->Write();
    matchedJPsiMomHist->Write();
    JPsiMomEff->Write();
    JPsiEtaMom->Write();
    matchedJPsiEtaMom->Write();
    matchedJPsiTrackDeltaR->Write();

    ofile->Close(); // Close output file
  }