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
#include <TGClient.h>
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

void kinematicPlots(){

    TString infile="eicReconOutput/EICreconOut_JPsiMuMu_1ifb_10x250ep_Pruned.root";
    std::string filename = infile.Data();
    std::string beam_config;
    double protonEnergy;
    std::string marker = "_Pruned.root";

    size_t pos = filename.rfind(marker);
    beam_config = filename.substr(pos - 8, 8);

    size_t x_pos = beam_config.find('x');
    std::string proton_part = beam_config.substr(x_pos + 1);  // e.g., "130ep"

    std::string proton_energy;
    for (char c : proton_part) {
        if (isdigit(c)) {
            proton_energy += c;
        } else {
            break;  // stop at first non-digit
        }
    }

    std::cout << "Proton energy: " << proton_energy << std::endl;

    protonEnergy = std::stod(proton_energy); 

    std::cout << "Extracted beam config: " << beam_config << std::endl;

    std::string outfilename = "outputs/kinematicsOutput_" + beam_config + ".root";

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
    TTreeReaderArray<double> partMass(tree_reader, "MCParticles.mass");
    TTreeReaderArray<unsigned int> partParb(tree_reader, "MCParticles.parents_begin");
    TTreeReaderArray<unsigned int> partPare(tree_reader, "MCParticles.parents_end");

    TTreeReaderArray<int> partParI(tree_reader, "_MCParticles_parents.index");

    // Get Reconstructed Track Information
    TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<float> trackEng(tree_reader, "ReconstructedChargedParticles.energy");
    TTreeReaderArray<int> trackPDG(tree_reader, "ReconstructedChargedParticles.PDG");
    TTreeReaderArray<float> trackMass(tree_reader, "ReconstructedChargedParticles.mass");
    TTreeReaderArray<float> trackCharge(tree_reader, "ReconstructedChargedParticles.charge");


    // Get Associations Between MCParticles and ReconstructedChargedParticles
    TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "ReconstructedChargedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations.simID");


    // Define Histograms
    TH1D* trueQ2 = new TH1D("trueQ2","True Q2 Distribution",100,0.,40.);
    TH1D* reconQ2 = new TH1D("reconQ2","Reconstructed Q2 Distribution",100,0.,40.);

    TH1D* truet = new TH1D("truet","True t Distribution",100,0.,2.);
    TH1D* reconT = new TH1D("reconT","Reconstructed t Distribution",100,0.,2.);

    TH1D* truexbjk = new TH1D("truexbjk","True xbj Distribution",100,0.,0.2);
    TH1D* reconXbjk = new TH1D("reconXbjk","Reconstructed xbj Distribution",100,0.,0.2);


    //Defining the four vectors
    ROOT::Math::PxPyPzEVector elec_beam; // initialized the 4 vector for electron beam
    ROOT::Math::PxPyPzEVector prot_beam; // initialized the 4 vector for proton beam
    ROOT::Math::PxPyPzEVector elec_mc; // initialized the 4 vector for truth electron
    ROOT::Math::PxPyPzEVector prot_mc; // initialized the 4 vector for truth electron
    ROOT::Math::PxPyPzEVector virtphoton_truth; // intialized the 4 vector for truth virtual photon
    ROOT::Math::PxPyPzEVector JPsi_mc; // initialized the 4 vector for truth J/PSi
    ROOT::Math::PxPyPzEVector muminus_mc; // initialized the 4 vector for truth muon minus
    ROOT::Math::PxPyPzEVector muplus_mc; // initialized the 4 vector for truth muon plus


    ROOT::Math::PxPyPzEVector elec_rec; // initialized the 4 vector for reconstructed electron
    ROOT::Math::PxPyPzEVector prot_rec; // initialized the 4 vector for reconstructed proton
    ROOT::Math::PxPyPzEVector virtphoton_rec; //intialized the 4 vector for reconstructed virtual photon
    ROOT::Math::PxPyPzEVector JPsi_rec; // initialized the 4 vector for reconstructed J/PSi
    ROOT::Math::PxPyPzEVector muminus_rec; // initialized the 4 vector for reconstructed muon minus
    ROOT::Math::PxPyPzEVector muplus_rec; // initialized the 4 vector for reconstructed muon plus

    double Q2_truth, Q2_rec;
    double t_truth, t_rec;
    double xbjk_truth, xbjk_rec;
    double partEng;

    bool scatE = false; // Flag to check if scattered electron is found

    // Defining initial colliding beams
    double eMass = 0.000510998950; //electron beam
    double eEng = 10;
    double e_pmag = sqrt(pow(eEng,2)-pow(eMass,2));
    double e_p1 = 0.;
    double e_p2 = 0.;
    double e_p3 = -1*e_pmag;
    elec_beam.SetPxPyPzE(e_p1, e_p2, e_p3, eEng); 
  
                
    double pMass = 0.93827208816; // proton beam
    double pEng = protonEnergy; //change
    double p_pmag = sqrt(pow(pEng,2)-pow(pMass,2));
    double c_a = 0.025;
    double p_p1 = -p_pmag*sin(c_a);
    double p_p2 = 0.;
    double p_p3 = p_pmag*cos(c_a);
    prot_beam.SetPxPyPzE(p_p1, p_p2, p_p3, pEng);


    while(tree_reader.Next()) { // Loop over events
      // Reset 4-vectors
      elec_mc.SetPxPyPzE(0.,0.,0.,0.);
      virtphoton_truth.SetPxPyPzE(0.,0.,0.,0.);
      JPsi_mc.SetPxPyPzE(0.,0.,0.,0.);
      muminus_mc.SetPxPyPzE(0.,0.,0.,0.);
      muplus_mc.SetPxPyPzE(0.,0.,0.,0.);
      elec_rec.SetPxPyPzE(0.,0.,0.,0.);
      virtphoton_rec.SetPxPyPzE(0.,0.,0.,0.);
      JPsi_rec.SetPxPyPzE(0.,0.,0.,0.);
      muminus_rec.SetPxPyPzE(0.,0.,0.,0.);
      muplus_rec.SetPxPyPzE(0.,0.,0.,0.);

      // Reset kinematic variables
      Q2_truth = 0.;
      Q2_rec = 0.;
      t_truth = 0.;
      t_rec = 0.;
      xbjk_truth = 0.;
      xbjk_rec = 0.;


      for(unsigned int i=0; i<partGenStat.GetSize(); i++) { // Loop over thrown particles
        partEng = sqrt(pow(partMomX[i],2) + pow(partMomY[i],2) + pow(partMomZ[i],2) + pow(partMass[i],2)); // Energy of all Monte Carlo particles
        
        if(partGenStat[i] == 1 && partPdg[i] == 11) { // Select stable thrown particles and look at electron
          elec_mc.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);
          for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
          {
            scatE = false;
            if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
            {
              for (int k = partParb[simuAssoc[j]]; k < partPare[simuAssoc[j]]; k++)
                {
                  int parentIndex = partParI[k];
                  if (partGenStat[parentIndex] == 4 && partPdg[parentIndex] == 11) // Check if the parent is a electron
                  {
                    scatE = true; // Set flag if a scattered electron is found
                  }
                  else
                  {
                    for (int l = partParb[parentIndex]; l < partPare[parentIndex]; l++)
                    {
                      int grandparentIndex = partParI[l];
                      if (partGenStat[grandparentIndex] == 4 && partPdg[grandparentIndex] == 11) // Check if the grandparent is a electron
                      {
                        scatE = true; // Set flag if a scattered electron is found
                      }
                    }
                  }
                }

              if (scatE)
              {
                elec_rec.SetPxPyPzE(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]], trackEng[recoAssoc[j]]);
                break;
              }
            }
          }
        }
        
        if(partGenStat[i] == 1 && partPdg[i] == 13) { // Look at mu minus
          muminus_mc.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);
          for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
          {
            if(simuAssoc[j] == i){
              muminus_rec.SetPxPyPzE(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]], trackEng[recoAssoc[j]]);
            }
          }
        }
          
        if(partGenStat[i] == 1 && partPdg[i] == -13) { // Look at mu plus
          muplus_mc.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);
          for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
          {
            if(simuAssoc[j] == i){
              muplus_rec.SetPxPyPzE(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]], trackEng[recoAssoc[j]]);
            }
          }
        }			 

        if(muminus_mc.P() > 0. && muplus_mc.P() > 0.) { // Check if both muons have momentum
          JPsi_mc = muminus_mc + muplus_mc; // Add the two muon 4-vectors to get the J/PSi 4-vector
          if (muminus_rec.P() > 0. && muplus_rec.P() > 0.) { // Check if both reconstructed muons have momentum
            JPsi_rec = muminus_rec + muplus_rec; // Add the two reconstructed muon 4-vectors to get the reconstructed J/PSi 4-vector
          }
        }

  
      } // for over thrown particles

      // Truth kinematic variables

      virtphoton_truth = (elec_beam - elec_mc);
      Q2_truth = -1*(virtphoton_truth.mag2());
      t_truth = -1*((virtphoton_truth - JPsi_mc).mag2());

      xbjk_truth = Q2_truth / (2* prot_beam.Dot(elec_beam - elec_mc)); 
      
      trueQ2->Fill(Q2_truth);
      truet->Fill(t_truth); 
      truexbjk->Fill(xbjk_truth);

      // Reconstructed kinematic variables

      if (elec_rec.P()>0. && muminus_rec.P() > 0. && muplus_rec.P() > 0.){ // if e', mu+, and mu- are in coincidence
  
        virtphoton_rec = (elec_beam - elec_rec);
        Q2_rec = -1*(virtphoton_rec.mag2()); 
        t_rec = -1*((virtphoton_rec - JPsi_rec).mag2());
        xbjk_rec = Q2_rec / (2* prot_beam.Dot(elec_beam - elec_rec)); 

        reconQ2->Fill(Q2_rec); 
        reconT->Fill(t_rec); 
        reconXbjk->Fill(xbjk_rec);
      }

    } 

    // Write histograms to file
    ofile->cd();
    ofile->mkdir("Kinematics");
    ofile->cd("Kinematics");
    trueQ2->Write();
    reconQ2->Write();
    truet->Write();
    reconT->Write();
    truexbjk->Write(); 
    reconXbjk->Write();
    ofile->cd("..");

    ofile->Close(); // Close output file
  }