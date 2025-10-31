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

void particleSearch()
{

    TString infile="reconOut/proton_20GeV_reconOut.root";

    double protonEnergy = 130;

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
    TTreeReaderArray<float> partCharge(tree_reader, "MCParticles.charge");
    TTreeReaderArray<unsigned int> partParb(tree_reader, "MCParticles.parents_begin");
    TTreeReaderArray<unsigned int> partPare(tree_reader, "MCParticles.parents_end");
    TTreeReaderArray<int> partParI(tree_reader, "_MCParticles_parents.index");

    // Get Reconstructed Track Information
    TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<int> trackPDG(tree_reader, "ReconstructedChargedParticles.PDG");
    TTreeReaderArray<float> trackMass(tree_reader, "ReconstructedChargedParticles.mass");
    TTreeReaderArray<float> trackCharge(tree_reader, "ReconstructedChargedParticles.charge");
    TTreeReaderArray<float> trackEng(tree_reader, "ReconstructedChargedParticles.energy");

    // Get Truth Seeded Reconstructed Track Information
    TTreeReaderArray<float> ttrackMomX(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x");
    TTreeReaderArray<float> ttrackMomY(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y");
    TTreeReaderArray<float> ttrackMomZ(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z");
    TTreeReaderArray<int> ttrackPDG(tree_reader, "ReconstructedTruthSeededChargedParticles.PDG");
    TTreeReaderArray<float> ttrackMass(tree_reader, "ReconstructedTruthSeededChargedParticles.mass");
    TTreeReaderArray<float> ttrackCharge(tree_reader, "ReconstructedTruthSeededChargedParticles.charge");
    TTreeReaderArray<float> ttrackEng(tree_reader, "ReconstructedTruthSeededChargedParticles.energy");

    // Get Associations Between MCParticles and ReconstructedChargedParticles
    TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "ReconstructedChargedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations.simID");

    // Get B0 Information
    TTreeReaderArray<unsigned int> recoAssocB0(tree_reader, "B0ECalClusterAssociations.recID");
    TTreeReaderArray<unsigned int> simuAssocB0(tree_reader, "B0ECalClusterAssociations.simID");
    TTreeReaderArray<float> B0Eng(tree_reader, "B0ECalClusters.energy");
    TTreeReaderArray<float> B0z(tree_reader, "B0ECalClusters.position.z");

    // Get Forward Detector Information
    TTreeReaderArray<float> RPEng(tree_reader, "ForwardRomanPotRecParticles.energy");
    TTreeReaderArray<int> RPpdg(tree_reader, "ForwardRomanPotRecParticles.PDG");
    TTreeReaderArray<float> RPMomX(tree_reader, "ForwardRomanPotRecParticles.momentum.x");
    TTreeReaderArray<float> RPMomY(tree_reader, "ForwardRomanPotRecParticles.momentum.y");
    TTreeReaderArray<float> RPMomZ(tree_reader, "ForwardRomanPotRecParticles.momentum.z");

    TTreeReaderArray<float> OffMEng(tree_reader, "ForwardOffMRecParticles.energy");

    // Ecal Information
    TTreeReaderArray<unsigned int> simuAssocEcalBarrel(tree_reader, "EcalBarrelClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocEcalBarrel(tree_reader, "EcalBarrelClusterAssociations.recID");
    TTreeReaderArray<float> EcalBarrelEng(tree_reader, "EcalBarrelClusters.energy");

    TTreeReaderArray<unsigned int> simuAssocEcalEndcapP(tree_reader, "EcalEndcapPClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocEcalEndcapP(tree_reader, "EcalEndcapPClusterAssociations.recID");    
    TTreeReaderArray<float> EcalEndcapPEng(tree_reader, "EcalEndcapPClusters.energy");

    TTreeReaderArray<unsigned int> simuAssocEcalEndcapN(tree_reader, "EcalEndcapNClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocEcalEndcapN(tree_reader, "EcalEndcapNClusterAssociations.recID");
    TTreeReaderArray<float> EcalEndcapNEng(tree_reader, "EcalEndcapNClusters.energy");

    // Hcal Information
    TTreeReaderArray<unsigned int> simuAssocHcalBarrel(tree_reader, "HcalBarrelClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocHcalBarrel(tree_reader, "HcalBarrelClusterAssociations.recID");
    TTreeReaderArray<float> HcalBarrelEng(tree_reader, "HcalBarrelClusters.energy");

    TTreeReaderArray<unsigned int> simuAssocHcalEndcapP(tree_reader, "HcalEndcapPInsertClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocHcalEndcapP(tree_reader, "HcalEndcapPInsertClusterAssociations.recID");    
    TTreeReaderArray<float> HcalEndcapPEng(tree_reader, "HcalEndcapPInsertClusters.energy");

    TTreeReaderArray<unsigned int> simuAssocHcalEndcapN(tree_reader, "HcalEndcapNClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocHcalEndcapN(tree_reader, "HcalEndcapNClusterAssociations.recID");
    TTreeReaderArray<float> HcalEndcapNEng(tree_reader, "HcalEndcapNClusters.energy");

    // Define 3 and 4-momentum vectors for particles

    TVector3 beamEMom;
    ROOT::Math::PxPyPzEVector beamE4Mom;
    TVector3 scatEMomT;
    ROOT::Math::PxPyPzEVector scatE4MomT;
    TVector3 scatEMomR;
    ROOT::Math::PxPyPzEVector scatE4MomR;

    TVector3 beampMom;
    ROOT::Math::PxPyPzEVector beamp4Mom;
    TVector3 scatpMomT;
    ROOT::Math::PxPyPzEVector scatp4MomT;
    TVector3 scatpMomR;
    ROOT::Math::PxPyPzEVector scatp4MomR;

    TVector3 muPlusMomT;
    ROOT::Math::PxPyPzEVector muPlus4MomT;
    TVector3 muPlusMomR;
    ROOT::Math::PxPyPzEVector muPlus4MomR;
    TVector3 muMinusMomT;
    ROOT::Math::PxPyPzEVector muMinus4MomT;
    TVector3 muMinusMomR;
    ROOT::Math::PxPyPzEVector muMinus4MomR;

    TVector3 JPsiMomT;
    ROOT::Math::PxPyPzEVector JPsi4MomT;
    TVector3 JPsiMomR;
    ROOT::Math::PxPyPzEVector JPsi4MomR;

    double Q2_truth, Q2_rec;
    double t_truth, t_rec;
    double xbjk_truth, xbjk_rec;
    double partEng;

    // Defining initial colliding beams
    double eMass = 0.000510998950; //electron beam
    double eEng = 10;
    double e_pmag = sqrt(pow(eEng,2)-pow(eMass,2));
    double e_p1 = 0.;
    double e_p2 = 0.;
    double e_p3 = -1*e_pmag;
    beamEMom = TVector3(e_p1, e_p2, e_p3);
    beamE4Mom  = ROOT::Math::PxPyPzEVector(e_p1, e_p2, e_p3, eEng); 
                
    double pMass = 0.93827208816; // proton beam
    double pEng = protonEnergy; //change
    double p_pmag = sqrt(pow(pEng,2)-pow(pMass,2));
    double c_a = 0.025;
    double p_p1 = -p_pmag*sin(c_a);
    double p_p2 = 0.;
    double p_p3 = p_pmag*cos(c_a);
    beampMom = TVector3(p_p1, p_p2, p_p3);
    beamp4Mom = ROOT::Math::PxPyPzEVector(p_p1, p_p2, p_p3, pEng);

    bool scatE = false; // Flag to check if scattered electron is found
    bool scatp = false; // Flag to check if scattered proton is found
    bool muPlus = false; // Flag to check if mu+ is found
    bool muMinus = false; // Flag to check if mu- is found
    int eventID = 0; // Event ID counter

    while(tree_reader.Next()) // Loop over events
    {
      eventID++;

      if (eventID % 10 == 0)
      {
        fprintf (stderr, "%4.2f Percent\r ", eventID*100.0/mychain->GetEntries());
        fflush (stderr);
      }

      // Reset vectors
      scatEMomT = TVector3(0.,0.,0.);
      scatE4MomT.SetPxPyPzE(0.,0.,0.,0.);
      scatEMomR = TVector3(0.,0.,0.);
      scatE4MomR.SetPxPyPzE(0.,0.,0.,0.);

      scatpMomT = TVector3(0.,0.,0.);
      scatp4MomT.SetPxPyPzE(0.,0.,0.,0.);
      scatpMomR = TVector3(0.,0.,0.);
      scatp4MomR.SetPxPyPzE(0.,0.,0.,0.);

      muPlusMomT = TVector3(0.,0.,0.);
      muPlus4MomT.SetPxPyPzE(0.,0.,0.,0.);
      muPlusMomR = TVector3(0.,0.,0.);
      muPlus4MomR.SetPxPyPzE(0.,0.,0.,0.);

      muMinusMomT = TVector3(0.,0.,0.);
      muMinus4MomT.SetPxPyPzE(0.,0.,0.,0.);
      muMinusMomR = TVector3(0.,0.,0.);
      muMinus4MomR.SetPxPyPzE(0.,0.,0.,0.);

      JPsiMomT = TVector3(0.,0.,0.);
      JPsi4MomT.SetPxPyPzE(0.,0.,0.,0.);
      JPsiMomR = TVector3(0.,0.,0.);
      JPsi4MomR.SetPxPyPzE(0.,0.,0.,0.);

      for(unsigned int i=0; i<partGenStat.GetSize(); i++)
      {       
        if (partGenStat[i] == 1)
        {
            partEng = sqrt(pow(partMomX[i],2) + pow(partMomY[i],2) + pow(partMomZ[i],2) + pow(partMass[i],2)); // Energy of Monte Carlo particle

            if (partPdg[i] == 2212)
            {
                scatpMomT = TVector3(partMomX[i],partMomY[i],partMomZ[i]);
                scatp4MomT.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);
            }
        }
      }
      
        for (int i = 0; i < trackPDG.GetSize(); i++)
        {
            if (i >= simuAssoc.GetSize() || simuAssoc[i] >= partPdg.GetSize()) continue; // Check to make sure the association index is valid           

            std::cout << "Track " << i  << "-- Simulated Associations: " << simuAssoc[i] << ", Reconstructed Association: " << recoAssoc[i] << std::endl;
            TVector3 recoMom(trackMomX[i],trackMomY[i],trackMomZ[i]);
            std::cout << "PDG = " << trackPDG[i] << ", eta = " << recoMom.PseudoRapidity() << ", E = " << trackEng[i] << ", m = " << std::sqrt(TMath::Abs(trackEng[i]*trackEng[i] - recoMom.Mag2())) << std::endl;

            for (int jB = 0; jB < EcalBarrelEng.GetSize(); jB++) // Look for associations in the Ecal Barrel
            {
                if (simuAssocEcalBarrel[jB] == simuAssoc[i])
                {
                    std::cout << "Ecal Barrel Cluster with energy: " << EcalBarrelEng[jB] << " GeV" << std::endl;
                }
            }
            for (int jP = 0; jP < EcalEndcapPEng.GetSize(); jP++) // Look for associations in the Ecal Endcap P
            {
                if (simuAssocEcalEndcapP[jP] == simuAssoc[i])
                {
                    std::cout << "Ecal Endcap P Cluster with energy: " << EcalEndcapPEng[jP] << " GeV" << std::endl;
                }
            }
            for (int jN = 0; jN < EcalEndcapNEng.GetSize(); jN++) // Look for associations in the Ecal Endcap N
            {
                if (simuAssocEcalEndcapN[jN] == simuAssoc[i])
                {
                    std::cout << "Ecal Endcap N Cluster with energy: " << EcalEndcapNEng[jN] << " GeV" << std::endl;
                }
            }

            if (HcalBarrelEng.GetSize() != simuAssocHcalBarrel.GetSize()) continue;
            for (int jB = 0; jB < simuAssocHcalBarrel.GetSize(); jB++) // Look for associations in the Ecal Barrel
            {
                if (simuAssocHcalBarrel[jB] == simuAssoc[i])
                {
                    std::cout << "Hcal Barrel Cluster with energy: " << HcalBarrelEng[jB] << " GeV" << std::endl;
                }
            }
            if (HcalEndcapPEng.GetSize() != simuAssocHcalEndcapP.GetSize()) continue;
            for (int jP = 0; jP < simuAssocHcalEndcapP.GetSize(); jP++) // Look for associations in the Ecal Endcap P
            {
                if (simuAssocHcalEndcapP[jP] == simuAssoc[i])
                {
                    std::cout << "Hcal Endcap P Cluster with energy: " << HcalEndcapPEng[jP] << " GeV" << std::endl;
                    exit(0);

                }
            }
            if (HcalEndcapNEng.GetSize() != simuAssocHcalEndcapN.GetSize()) continue;
            for (int jN = 0; jN < simuAssocHcalEndcapN.GetSize(); jN++) // Look for associations in the Ecal Endcap N
            {
                if (simuAssocHcalEndcapN[jN] == simuAssoc[i])
                {
                    std::cout << "Hcal Endcap N Cluster with energy: " << HcalEndcapNEng[jN] << " GeV" << std::endl;
                }
            }

            TVector3 trueMom(partMomX[simuAssoc[i]],partMomY[simuAssoc[i]],partMomZ[simuAssoc[i]]);
            std::cout << "True PDG = " << partPdg[simuAssoc[i]] << ", eta = " << trueMom.PseudoRapidity() << ", E = " << sqrt(trueMom.Mag2() + pow(partMass[simuAssoc[i]],2)) << std::endl;

            std::cout << " " << std::endl;

        }

    }


    std::cout << "" << std::endl;
    std::cout << "Event Processing Complete" << std::endl;
}