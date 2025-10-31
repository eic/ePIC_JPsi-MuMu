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

void calPlots()
{

    TString infileMuMinus="reconOut/muMinus_20GeV_reconOut.root";
    TString infileMuPlus="reconOut/muPlus_20GeV_reconOut.root";
    std::string outfilename = "outputs/mumuCalPlotsOutput.root";

    // Set output file for the histograms
    TFile *ofile = TFile::Open(outfilename.c_str(),"RECREATE");

    // Set up input file chain
    TChain *mychain = new TChain("events");
    mychain->Add(infileMuMinus);
    //mychain->Add(infileMuPlus);

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

    TTreeReaderArray<float> EcalBarrelHitsEng(tree_reader, "EcalBarrelImagingRecHits.energy");
    TTreeReaderArray<float> EcalBarrelHitsX(tree_reader, "EcalBarrelImagingRecHits.position.x");
    TTreeReaderArray<float> EcalBarrelHitsY(tree_reader, "EcalBarrelImagingRecHits.position.y");
    TTreeReaderArray<float> EcalBarrelHitsZ(tree_reader, "EcalBarrelImagingRecHits.position.z");


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

    TTreeReaderArray<float> HcalBarrelHitsEng(tree_reader, "HcalBarrelRecHits.energy");
    TTreeReaderArray<float> HcalBarrelHitsX(tree_reader, "HcalBarrelRecHits.position.x");
    TTreeReaderArray<float> HcalBarrelHitsY(tree_reader, "HcalBarrelRecHits.position.y");
    TTreeReaderArray<float> HcalBarrelHitsZ(tree_reader, "HcalBarrelRecHits.position.z");


    // Define histograms

    TH2D *hEcalBarrelHitsXY = new TH2D("hEcalBarrelHitsXY","Ecal Barrel Hits XY; x (cm); y (cm)",125,-2500.,2500.,125,-2500.,2500.);
    TH2D *hEcalBarrelHitsZY = new TH2D("hEcalBarrelHitsZY","Ecal Barrel Hits ZY; z (cm); y (cm)",150,-3000.,3000.,125,-2500.,2500.);
    TH2D *hEcalBarrelHitsZX = new TH2D("hEcalBarrelHitsZX","Ecal Barrel Hits ZX; z (cm); x (cm)",150,-3000.,3000.,125,-2500.,2500.);

    TH2D *hEcalBarrelHitsXYevent = new TH2D("hEcalBarrelHitsXYevent","Ecal Barrel Hits XY; x (cm); y (cm)",125,-2500.,2500.,125,-2500.,2500.);
    TH2D *hEcalBarrelHitsZYevent = new TH2D("hEcalBarrelHitsZYevent","Ecal Barrel Hits ZY; z (cm); y (cm)",150,-3000.,3000.,125,-2500.,2500.);
    TH2D *hEcalBarrelHitsZXevent = new TH2D("hEcalBarrelHitsZXevent","Ecal Barrel Hits ZX; z (cm); x (cm)",150,-3000.,3000.,125,-2500.,2500.);

    TH2D *hHcalBarrelHitsXY = new TH2D("hHcalBarrelHitsXY","Hcal Barrel Hits XY; x (cm); y (cm)",125,-2500.,2500.,125,-2500.,2500.);
    TH2D *hHcalBarrelHitsZY = new TH2D("hHcalBarrelHitsZY","Hcal Barrel Hits ZY; z (cm); y (cm)",150,-3000.,3000.,125,-2500.,2500.);
    TH2D *hHcalBarrelHitsZX = new TH2D("hHcalBarrelHitsZX","Hcal Barrel Hits ZX; z (cm); x (cm)",150,-3000.,3000.,125,-2500.,2500.);

    TH2D *hHcalBarrelHitsXYevent = new TH2D("hHcalBarrelHitsXYevent","Hcal Barrel Hits XY; x (cm); y (cm)",125,-2500.,2500.,125,-2500.,2500.);
    TH2D *hHcalBarrelHitsZYevent = new TH2D("hHcalBarrelHitsZYevent","Hcal Barrel Hits ZY; z (cm); y (cm)",150,-3000.,3000.,125,-2500.,2500.);
    TH2D *hHcalBarrelHitsZXevent = new TH2D("hHcalBarrelHitsZXevent","Hcal Barrel Hits ZX; z (cm); x (cm)",150,-3000.,3000.,125,-2500.,2500.);

    int eventID = 0;

    while(tree_reader.Next()) // Loop over events
    {
      eventID++;

      hEcalBarrelHitsXYevent->Reset();
      hEcalBarrelHitsZYevent->Reset();
      hEcalBarrelHitsZXevent->Reset(); 
      hHcalBarrelHitsXYevent->Reset();
      hHcalBarrelHitsZYevent->Reset();
      hHcalBarrelHitsZXevent->Reset();

      if (eventID % 10 == 0)
      {
        fprintf (stderr, "%4.2f Percent\r ", eventID*100.0/mychain->GetEntries());
        fflush (stderr);
      }

      for (int iEH = 0; iEH < EcalBarrelHitsEng.GetSize(); iEH++)
      {
        hEcalBarrelHitsXY->Fill(EcalBarrelHitsX[iEH],EcalBarrelHitsY[iEH],EcalBarrelHitsEng[iEH]);
        hEcalBarrelHitsZY->Fill(EcalBarrelHitsZ[iEH],EcalBarrelHitsY[iEH],EcalBarrelHitsEng[iEH]);
        hEcalBarrelHitsZX->Fill(EcalBarrelHitsZ[iEH],EcalBarrelHitsX[iEH],EcalBarrelHitsEng[iEH]);
        hEcalBarrelHitsXYevent->Fill(EcalBarrelHitsX[iEH],EcalBarrelHitsY[iEH],EcalBarrelHitsEng[iEH]);
        hEcalBarrelHitsZYevent->Fill(EcalBarrelHitsZ[iEH],EcalBarrelHitsY[iEH],EcalBarrelHitsEng[iEH]);
        hEcalBarrelHitsZXevent->Fill(EcalBarrelHitsZ[iEH],EcalBarrelHitsX[iEH],EcalBarrelHitsEng[iEH]);
      }

        for (int iHH = 0; iHH < HcalBarrelHitsEng.GetSize(); iHH++)
        {
          hHcalBarrelHitsXY->Fill(HcalBarrelHitsX[iHH],HcalBarrelHitsY[iHH],HcalBarrelHitsEng[iHH]);
          hHcalBarrelHitsZY->Fill(HcalBarrelHitsZ[iHH],HcalBarrelHitsY[iHH],HcalBarrelHitsEng[iHH]);
          hHcalBarrelHitsZX->Fill(HcalBarrelHitsZ[iHH],HcalBarrelHitsX[iHH],HcalBarrelHitsEng[iHH]);
          hHcalBarrelHitsXYevent->Fill(HcalBarrelHitsX[iHH],HcalBarrelHitsY[iHH],HcalBarrelHitsEng[iHH]);
          hHcalBarrelHitsZYevent->Fill(HcalBarrelHitsZ[iHH],HcalBarrelHitsY[iHH],HcalBarrelHitsEng[iHH]);
          hHcalBarrelHitsZXevent->Fill(HcalBarrelHitsZ[iHH],HcalBarrelHitsX[iHH],HcalBarrelHitsEng[iHH]);
        }

        if(EcalBarrelHitsEng.GetSize() > 10 && HcalBarrelHitsEng.GetSize() > 10)
        {
          TCanvas *c1 = new TCanvas("c1","Ecal and Hcal Barrel Hits",1200,800);
          c1->Divide(2,3);
          c1->cd(1);
          hEcalBarrelHitsXYevent->Draw("COLZ");
          c1->cd(2);
          hHcalBarrelHitsXYevent->Draw("COLZ");
          c1->cd(3);
          hEcalBarrelHitsZXevent->Draw("COLZ");
          c1->cd(4);
          hHcalBarrelHitsZXevent->Draw("COLZ");
          c1->cd(5);
          hEcalBarrelHitsZYevent->Draw("COLZ");
          c1->cd(6);
          hHcalBarrelHitsZYevent->Draw("COLZ");
          c1->Update();
          std::cout << "Press Enter to continue..." << std::endl;
          std::cin.get();
          delete c1;
        }

    }

    // Write output histograms to file
    ofile->cd();
    hEcalBarrelHitsXY->Write();
    hEcalBarrelHitsZY->Write();
    hEcalBarrelHitsZX->Write();
    hHcalBarrelHitsXY->Write();
    hHcalBarrelHitsZY->Write();
    hHcalBarrelHitsZX->Write();
    ofile->Close();
}