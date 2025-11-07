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
#include "TPie.h"

void epPlots()
{

    //TString infile="../eicReconOutput/EICreconOut_JPsiMuMu_10ifb_10x130ep_Pruned.root";
    TString infile="reconOut/mu-pi-ka-e-p_20GeV_reconOut.root";

    double protonEnergy = 130; 

    //std::string outfilename = "outputs/JPsiMuMu_10ifb_10x130ep_epPlots.root";
    std::string outfilename = "outputs/mu-pi-ka-e-p_epPlots.root";

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
    TTreeReaderArray<float> partCharge(tree_reader, "MCParticles.charge");
    TTreeReaderArray<unsigned int> partParb(tree_reader, "MCParticles.parents_begin");
    TTreeReaderArray<unsigned int> partPare(tree_reader, "MCParticles.parents_end");
    TTreeReaderArray<int> partParI(tree_reader, "_MCParticles_parents.index");

    // Get Reconstructed Track Information
    TTreeReaderArray<float> recoMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> recoMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> recoMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<int> trackPDG(tree_reader, "ReconstructedChargedParticles.PDG");
    TTreeReaderArray<float> trackMass(tree_reader, "ReconstructedChargedParticles.mass");
    TTreeReaderArray<float> trackCharge(tree_reader, "ReconstructedChargedParticles.charge");
    TTreeReaderArray<float> trackEng(tree_reader, "ReconstructedChargedParticles.energy");

    // Get Truth Seeded Reconstructed Track Information
    TTreeReaderArray<float> trecoMomX(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x");
    TTreeReaderArray<float> trecoMomY(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y");
    TTreeReaderArray<float> trecoMomZ(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z");
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

    // Pseudo-rapidity plots by particle
    TH1D *protonRecEta = new TH1D("protonEta","Proton Pseudo-rapidity; eta; Counts",100,-5.,5.);
    TH1D *electronRecEta = new TH1D("electronEta","Electron Pseudo-rapidity; eta; Counts",100,-5.,5.);
    TH1D *pionRecEta = new TH1D("pionEta","Pion Pseudo-rapidity; eta; Counts",100,-5.,5.);
    TH1D *kaonRecEta = new TH1D("kaonEta","Kaon Pseudo-rapidity; eta; Counts",100,-5.,5.);
    TH1D *muonRecEta = new TH1D("muonEta","Muon Pseudo-rapidity; eta; Counts",100,-5.,5.);

    // Momentum plots by particle
    TH1D *protonRecP = new TH1D("protonP","Proton Momentum; p [GeV/c]; Counts",100,0.,25.);
    TH1D *electronRecP = new TH1D("electronP","Electron Momentum; p [GeV/c]; Counts",100,0.,25.);
    TH1D *pionRecP = new TH1D("pionP","Pion Momentum; p [GeV/c]; Counts",100,0.,25.);
    TH1D *kaonRecP = new TH1D("kaonP","Kaon Momentum; p [GeV/c]; Counts",100,0.,25.);
    TH1D *muonRecP = new TH1D("muonP","Muon Momentum; p [GeV/c]; Counts",100,0.,25.);

    // Ecal ep Plots by particle
    TH2D *protonEpTrueEcal = new TH2D("protonEpTrueEcal","Ecal Energy vs Track Momentum for Protons; p;  E/p",100,0.,25.,200,0.,2.0);
    TH2D *electronEpTrueEcal = new TH2D("electronEpTrueEcal","Ecal Energy vs Track Momentum for Electrons; p;  E/p",100,0.,25.,100,0.,2.0);
    TH2D *pionEpTrueEcal = new TH2D("pionEpTrueEcal","Ecal Energy vs Track Momentum for Pions; p;  E/p",100,0.,25.,100,0.,2.0);
    TH2D *kaonEpTrueEcal = new TH2D("kaonEpTrueEcal","Ecal Energy vs Track Momentum for Kaons; p;  E/p",100,0.,25.,100,0.,2.0);
    TH2D *muonEpTrueEcal = new TH2D("muonEpTrueEcal","Ecal Energy vs Track Momentum for Muons; p;  E/p",100,0.,25.,100,0.,2.0);

    TH2D *protonEpRecEcal = new TH2D("protonEpRecEcal","Ecal Energy vs Track Momentum for Protons; p;  E/p",100,0.,25.,200,0.,2.0);
    TH2D *electronEpRecEcal = new TH2D("electronEpRecEcal","Ecal Energy vs Track Momentum for Electrons; p;  E/p",100,0.,25.,100,0.,2.0);
    TH2D *pionEpRecEcal = new TH2D("pionEpRecEcal","Ecal Energy vs Track Momentum for Pions; p;  E/p",100,0.,25.,100,0.,2.0);
    TH2D *kaonEpRecEcal = new TH2D("kaonEpRecEcal","Ecal Energy vs Track Momentum for Kaons; p;  E/p",100,0.,25.,100,0.,2.0);
    TH2D *muonEpRecEcal = new TH2D("muonEpRecEcal","Ecal Energy vs Track Momentum for Muons; p;  E/p",100,0.,25.,100,0.,2.0);

    // Hcal ep Plots by particle
    TH2D *protonEpTrueHcal = new TH2D("protonEpTrueHcal","Hcal Energy vs Track Momentum for Protons; p;  E/p",100,0.,25.,200,0.,2.0);
    TH2D *electronEpTrueHcal = new TH2D("electronEpTrueHcal","Hcal Energy vs Track Momentum for Electrons; p;  E/p",100,0.,25.,100,0.,2.0);
    TH2D *pionEpTrueHcal = new TH2D("pionEpTrueHcal","Hcal Energy vs Track Momentum for Pions; p;  E/p",100,0.,25.,100,0.,2.0);
    TH2D *kaonEpTrueHcal = new TH2D("kaonEpTrueHcal","Hcal Energy vs Track Momentum for Kaons; p;  E/p",100,0.,25.,100,0.,2.0);
    TH2D *muonEpTrueHcal = new TH2D("muonEpTrueHcal","Hcal Energy vs Track Momentum for Muons; p;  E/p",100,0.,25.,100,0.,2.0);

    TH2D *protonEpRecHcal = new TH2D("protonEpRecHcal","Hcal Energy vs Track Momentum for Protons; p;  E/p",100,0.,25.,200,0.,2.0);
    TH2D *electronEpRecHcal = new TH2D("electronEpRecHcal","Hcal Energy vs Track Momentum for Electrons; p;  E/p",100,0.,25.,100,0.,2.0);
    TH2D *pionEpRecHcal = new TH2D("pionEpRecHcal","Hcal Energy vs Track Momentum for Pions; p;  E/p",100,0.,25.,100,0.,2.0);
    TH2D *kaonEpRecHcal = new TH2D("kaonEpRecHcal","Hcal Energy vs Track Momentum for Kaons; p;  E/p",100,0.,25.,100,0.,2.0);
    TH2D *muonEpRecHcal = new TH2D("muonEpRecHcal","Hcal Energy vs Track Momentum for Muons; p;  E/p",100,0.,25.,100,0.,2.0);

    // Ecal vs Hcal energy Plots by particle
    TH2D *protonEpTrueEHcal = new TH2D("protonEpTrueEHcal","Ecal Energy vs Hcal Energy for Protons; E/p;  E/p",100,0.,2.0,100,0.,2.0);
    TH2D *electronEpTrueEHcal = new TH2D("electronEpTrueEHcal","Ecal Energy vs Hcal Energy for Electrons; E/p;  E/p",100,0.,2.0,100,0.,2.0);
    TH2D *pionEpTrueEHcal = new TH2D("pionEpTrueEHcal","Ecal Energy vs Hcal Energy for Pions; E/p;  E/p",100,0.,2.0,100,0.,2.0);
    TH2D *kaonEpTrueEHcal = new TH2D("kaonEpTrueEHcal","Ecal Energy vs Hcal Energy for Kaons; E/p;  E/p",100,0.,2.0,100,0.,2.0);
    TH2D *muonEpTrueEHcal = new TH2D("muonEpTrueEHcal","Ecal Energy vs Hcal Energy for Muons; E/p;  E/p",100,0.,2.0,100,0.,2.0);

    TH2D *protonEpRecEHcal = new TH2D("protonEpRecEHcal","Ecal Energy vs Hcal Energy for Protons; E/p;  E/p",100,0.,2.0,100,0.,2.0);
    TH2D *electronEpRecEHcal = new TH2D("electronEpRecEHcal","Ecal Energy vs Hcal Energy for Electrons; E/p;  E/p",100,0.,2.0,100,0.,2.0);
    TH2D *pionEpRecEHcal = new TH2D("pionEpRecEHcal","Ecal Energy vs Hcal Energy for Pions; E/p;  E/p",100,0.,2.0,100,0.,2.0);
    TH2D *kaonEpRecEHcal = new TH2D("kaonEpRecEHcal","Ecal Energy vs Hcal Energy for Kaons; E/p;  E/p",100,0.,2.0,100,0.,2.0);
    TH2D *muonEpRecEHcal = new TH2D("muonEpRecEHcal","Ecal Energy vs Hcal Energy for Muons; E/p;  E/p",100,0.,2.0,100,0.,2.0);

    // Cluster Size Plots by particle
    TH1D *protonEcalClusterSize = new TH1D("protonEcalClusterSize","Ecal Barrel Cluster Size for Protons; Size; Counts",100,0.,50.);
    TH1D *electronEcalClusterSize = new TH1D("electronEcalClusterSize","Ecal Barrel Cluster Size for Electrons; Size; Counts",100,0.,50.);
    TH1D *pionEcalClusterSize = new TH1D("pionEcalClusterSize","Ecal Barrel Cluster Size for Pions; Size; Counts",100,0.,50.);
    TH1D *kaonEcalClusterSize = new TH1D("kaonEcalClusterSize","Ecal Barrel Cluster Size for Kaons; Size; Counts",100,0.,50.);
    TH1D *muonEcalClusterSize = new TH1D("muonEcalClusterSize","Ecal Barrel Cluster Size for Muons; Size; Counts",100,0.,50.);

    TH1D *protonHcalClusterSize = new TH1D("protonHcalClusterSize","Hcal Barrel Cluster Size for Protons; Size; Counts",100,0.,50.);
    TH1D *electronHcalClusterSize = new TH1D("electronHcalClusterSize","Hcal Barrel Cluster Size for Electrons; Size; Counts",100,0.,50.);
    TH1D *pionHcalClusterSize = new TH1D("pionHcalClusterSize","Hcal Barrel Cluster Size for Pions; Size; Counts",100,0.,50.);
    TH1D *kaonHcalClusterSize = new TH1D("kaonHcalClusterSize","Hcal Barrel Cluster Size for Kaons; Size; Counts",100,0.,50.);
    TH1D *muonHcalClusterSize = new TH1D("muonHcalClusterSize","Hcal Barrel Cluster Size for Muons; Size; Counts",100,0.,50.);

    TH1D *protonEHcalClusterSize = new TH1D("protonEHcalClusterSize","Ecal + Hcal Barrel Cluster Size for Protons; Size; Counts",200,0.,100.);
    TH1D *electronEHcalClusterSize = new TH1D("electronEHcalClusterSize","Ecal + Hcal Barrel Cluster Size for Electrons; Size; Counts",200,0.,100.);
    TH1D *pionEHcalClusterSize = new TH1D("pionEHcalClusterSize","Ecal + Hcal Barrel Cluster Size for Pions; Size; Counts",200,0.,100.);
    TH1D *kaonEHcalClusterSize = new TH1D("kaonEHcalClusterSize","Ecal + Hcal Barrel Cluster Size for Kaons; Size; Counts",200,0.,100.);
    TH1D *muonEHcalClusterSize = new TH1D("muonEHcalClusterSize","Ecal + Hcal Barrel Cluster Size for Muons; Size; Counts",200,0.,100.);

    // Post cut plot

    TH1D *initialPID = new TH1D("initialPID","Initial Particle ID Distribution; PDG; Counts",2500,0,2500);
    TH1D *trackPIDCutPID = new TH1D("trackPIDCutPID","Particle ID Distribution after Track PID Cut; PDG; Counts",2500,0,2500);
    TH1D *energyCutPID = new TH1D("energyCutPID","Particle ID Distribution after Energy and Mass Cut; PDG; Counts",2500,0,2500);
    TH1D *epCutPID = new TH1D("epCutPID","Particle ID Distribution after E/p Cut; PDG; Counts",2500,0,2500);
    TH1D *clusterSizePID = new TH1D("clusterSizePID","Particle ID Distribution after Cluster Size Cut; PDG; Counts",2500,0,2500);

    TH1D *pionPostCutEta = new TH1D("pionPCEta","Pion Pseudo-rapidity post cuts; eta; Counts",100,-5.,5.);
    TH1D *muonPostCutEta = new TH1D("muonPCEta","Muon Pseudo-rapidity post cuts; eta; Counts",100,-5.,5.);

    TH1D *pionPostCutP = new TH1D("pionPCP","Pion Momentum post cuts; p [GeV/c]; Counts",100,0.,25.);
    TH1D *muonPostCutP = new TH1D("muonPCP","Muon Momentum post cuts; p [GeV/c]; Counts",100,0.,25.);

    // Purity and Efficiency Counters

    double momentumBins[20] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5};
    double muonTrueCounts[20] = {0};
    double muonRecCounts[20] = {0};
    double totalParticleTrueCounts[20] = {0};
    double totalParticleRecCounts[20] = {0};
    double muonPurityTrue[20] = {0};
    double muonPurityRec[20] = {0};
    double muonEfficiency[20] = {0};

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

    int eventID = 0; // Event ID counter
    double ECalEnergy = 0.;
    double HCalEnergy = 0.;
    double EcalHits = 0.;
    double HcalHits = 0.;

    double trueParticles[5] = {0,0,0,0,0};
    double totalParticles[5] = {0,0,0,0,0};
    double muons[5] = {0,0,0,0,0};
    double pions[5] = {0,0,0,0,0};
    double kaons[5] = {0,0,0,0,0};
    double protons[5] = {0,0,0,0,0};
    double electrons[5] = {0,0,0,0,0};

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
            double partEta = TVector3(partMomX[i],partMomY[i],partMomZ[i]).Eta();

            if (partPdg[i]== 2212){
                scatpMomT = TVector3(partMomX[i],partMomY[i],partMomZ[i]);
                scatp4MomT.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);
            }

            if (partEta < -4 || partEta > 4) continue;
            if (TVector3(partMomX[i],partMomY[i],partMomZ[i]).Mag() < 0.05) continue;
                
            totalParticleTrueCounts[int((TVector3(partMomX[i],partMomY[i],partMomZ[i]).Mag()-0.5))] += 1;

            switch (TMath::Abs(partPdg[i]))
            {
                case 13:
                    trueParticles[0] += 1.;
                    muonTrueCounts[int((TVector3(partMomX[i],partMomY[i],partMomZ[i]).Mag()-0.5))] += 1;
                    break;
                case 211:
                    trueParticles[1] += 1.;
                    break;
                case 321:
                    trueParticles[2] += 1.;
                    break;
                case 2212:
                    trueParticles[3] += 1.;
                    break;
                case 11:
                    trueParticles[4] += 1.;
                    break;
                default:
                    break;
            }
        }
      }
      
        for (int i = 0; i < trackPDG.GetSize(); i++)
        {
            ECalEnergy = 0.;
            HCalEnergy = 0.;
            EcalHits = 0.;
            HcalHits = 0.;
            TVector3 recoMom(recoMomX[i],recoMomY[i],recoMomZ[i]);
            ROOT::Math::PxPyPzEVector reco4Mom(recoMomX[i],recoMomY[i],recoMomZ[i], trackEng[i]);

            if (recoMom.Eta() < -4 || recoMom.Eta() > 4) continue;
            if (recoMom.Mag() < 0.05) continue;

            if (i >= simuAssoc.GetSize() || simuAssoc[i] >= partPdg.GetSize()) continue;

            totalParticles[0] += 1.;

            if (EcalBarrelEng.GetSize() == simuAssocEcalBarrel.GetSize())
            {
                for (int jB = 0; jB < simuAssocEcalBarrel.GetSize(); jB++) // Look for associations in the Ecal Barrel
                {
                    if (simuAssocEcalBarrel[jB] == simuAssoc[i])
                    {
                        ECalEnergy += EcalBarrelEng[jB];
                        EcalHits += 1.;
                    }
                }
            }
            if (EcalEndcapPEng.GetSize() == simuAssocEcalEndcapP.GetSize()) 
            {
                for (int jP = 0; jP < simuAssocEcalEndcapP.GetSize(); jP++) // Look for associations in the Ecal Endcap P
                {
                    if (simuAssocEcalEndcapP[jP] == simuAssoc[i])
                    {
                        ECalEnergy += EcalEndcapPEng[jP];
                        EcalHits += 1.;
                    }
                }
            }
            if (EcalEndcapNEng.GetSize() == simuAssocEcalEndcapN.GetSize())
            {
                for (int jN = 0; jN < simuAssocEcalEndcapN.GetSize(); jN++) // Look for associations in the Ecal Endcap N
                {
                    if (simuAssocEcalEndcapN[jN] == simuAssoc[i])
                    {
                        ECalEnergy += EcalEndcapNEng[jN];
                        EcalHits += 1.;
                    }
                }
            }

            if (HcalBarrelEng.GetSize() == simuAssocHcalBarrel.GetSize())
            {
                for (int jB = 0; jB < simuAssocHcalBarrel.GetSize(); jB++) // Look for associations in the Ecal Barrel
                {
                    if (simuAssocHcalBarrel[jB] == simuAssoc[i])
                    {
                        HCalEnergy += HcalBarrelEng[jB];
                        HcalHits += 1.;
                    }
                }
            }
            
            if (HcalEndcapPEng.GetSize() == simuAssocHcalEndcapP.GetSize())
            {
                for (int jP = 0; jP < simuAssocHcalEndcapP.GetSize(); jP++) // Look for associations in the Ecal Endcap P
                {
                    if (simuAssocHcalEndcapP[jP] == simuAssoc[i])
                    {
                        HCalEnergy += HcalEndcapPEng[jP];
                        HcalHits += 1.;
                    }
                }
            }
            
            if (HcalEndcapNEng.GetSize() == simuAssocHcalEndcapN.GetSize())
            {
                for (int jN = 0; jN < simuAssocHcalEndcapN.GetSize(); jN++) // Look for associations in the Ecal Endcap N
                {
                    if (simuAssocHcalEndcapN[jN] == simuAssoc[i])
                    {
                        HCalEnergy += HcalEndcapNEng[jN];
                        HcalHits += 1.;
                    }
                }
            }

            
            TVector3 trueMom(partMomX[simuAssoc[i]],partMomY[simuAssoc[i]],partMomZ[simuAssoc[i]]);
            ROOT::Math::PxPyPzEVector true4Mom(partMomX[simuAssoc[i]],partMomY[simuAssoc[i]],partMomZ[simuAssoc[i]], partEng);
            double trueMass = std::sqrt( trueMom.Mag2() - partEng*partEng);

            //if (HCalEnergy/trueMom.Mag() > 2.0) std::cout << "HCal Ratio: " <<  HCalEnergy/trueMom.Mag() << std::endl;

            if (ECalEnergy > 0. && HCalEnergy > 0.)
            {
                switch (TMath::Abs(partPdg[simuAssoc[i]]))
                {
                    case 2212:
                        protonRecEta->Fill(recoMom.PseudoRapidity());
                        protonRecP->Fill(recoMom.Mag());
                        protonEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        protonEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        protonEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        protonEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        protonEpTrueEHcal->Fill(ECalEnergy/trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        protonEpRecEHcal->Fill(ECalEnergy/recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        protonEcalClusterSize->Fill(EcalHits);
                        protonHcalClusterSize->Fill(HcalHits);
                        protonEHcalClusterSize->Fill(EcalHits+HcalHits);
                        break;
                    case 11:
                        electronRecEta->Fill(recoMom.PseudoRapidity());
                        electronRecP->Fill(recoMom.Mag());
                        electronEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        electronEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        electronEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        electronEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        electronEpTrueEHcal->Fill(ECalEnergy/trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        electronEpRecEHcal->Fill(ECalEnergy/recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        electronEcalClusterSize->Fill(EcalHits);
                        electronHcalClusterSize->Fill(HcalHits);
                        electronEHcalClusterSize->Fill(EcalHits+HcalHits);
                        break;
                    case 211:
                        pionRecEta->Fill(recoMom.PseudoRapidity());
                        pionRecP->Fill(recoMom.Mag());
                        pionEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        pionEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        pionEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        pionEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        pionEpTrueEHcal->Fill(ECalEnergy/trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        pionEpRecEHcal->Fill(ECalEnergy/recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        pionEcalClusterSize->Fill(EcalHits);
                        pionHcalClusterSize->Fill(HcalHits);
                        pionEHcalClusterSize->Fill(EcalHits+HcalHits);
                        break;
                    case 321:
                        kaonRecEta->Fill(recoMom.PseudoRapidity());
                        kaonRecP->Fill(recoMom.Mag());
                        kaonEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        kaonEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        kaonEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        kaonEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        kaonEpTrueEHcal->Fill(ECalEnergy/trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        kaonEpRecEHcal->Fill(ECalEnergy/recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        kaonEcalClusterSize->Fill(EcalHits);
                        kaonHcalClusterSize->Fill(HcalHits);
                        kaonEHcalClusterSize->Fill(EcalHits+HcalHits);
                        break;
                    case 13:
                        muonRecEta->Fill(recoMom.PseudoRapidity());
                        muonRecP->Fill(recoMom.Mag());
                        muonEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        muonEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        muonEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        muonEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        muonEpTrueEHcal->Fill(ECalEnergy/trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        muonEpRecEHcal->Fill(ECalEnergy/recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        muonEcalClusterSize->Fill(EcalHits);
                        muonHcalClusterSize->Fill(HcalHits);
                        muonEHcalClusterSize->Fill(EcalHits+HcalHits);
                        break;
                    default:
                        break;
                }
            }
            else if (ECalEnergy > 0. && HCalEnergy == 0.)
            {
                switch (TMath::Abs(partPdg[simuAssoc[i]]))
                {
                    case 2212:
                        protonRecEta->Fill(recoMom.PseudoRapidity());
                        protonRecP->Fill(recoMom.Mag());
                        protonEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        protonEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        protonEcalClusterSize->Fill(EcalHits);
                        break;
                    case 11:
                        electronRecEta->Fill(recoMom.PseudoRapidity());
                        electronRecP->Fill(recoMom.Mag());
                        electronEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        electronEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        electronEcalClusterSize->Fill(EcalHits);
                        break;
                    case 211:
                        pionRecEta->Fill(recoMom.PseudoRapidity());
                        pionRecP->Fill(recoMom.Mag());
                        pionEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        pionEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        pionEcalClusterSize->Fill(EcalHits);
                        break;
                    case 321:
                        kaonRecEta->Fill(recoMom.PseudoRapidity());
                        kaonRecP->Fill(recoMom.Mag());
                        kaonEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        kaonEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        kaonEcalClusterSize->Fill(EcalHits);
                        break;
                    case 13:
                        muonRecEta->Fill(recoMom.PseudoRapidity());
                        muonRecP->Fill(recoMom.Mag());
                        muonEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        muonEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        muonEcalClusterSize->Fill(EcalHits);
                        break;
                    default:
                        break;
                }
            }
            else if (ECalEnergy == 0. && HCalEnergy > 0.)
            {
                switch (TMath::Abs(partPdg[simuAssoc[i]]))
                {
                    case 2212:
                        protonRecEta->Fill(recoMom.PseudoRapidity());
                        protonRecP->Fill(recoMom.Mag());
                        protonEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        protonEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        protonHcalClusterSize->Fill(HcalHits);
                        break;
                    case 11:
                        electronRecEta->Fill(recoMom.PseudoRapidity());
                        electronRecP->Fill(recoMom.Mag());
                        electronEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        electronEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        electronHcalClusterSize->Fill(HcalHits);
                        break;
                    case 211:
                        pionRecEta->Fill(recoMom.PseudoRapidity());
                        pionRecP->Fill(recoMom.Mag());
                        pionEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        pionEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        pionHcalClusterSize->Fill(HcalHits);
                        break;
                    case 321:
                        kaonRecEta->Fill(recoMom.PseudoRapidity());
                        kaonRecP->Fill(recoMom.Mag());
                        kaonEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        kaonEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        kaonHcalClusterSize->Fill(HcalHits);
                        break;
                    case 13:
                        muonRecEta->Fill(recoMom.PseudoRapidity());
                        muonRecP->Fill(recoMom.Mag());
                        muonEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        muonEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        muonHcalClusterSize->Fill(HcalHits);
                        break;
                    default:
                        break;
                }
            }
            else
            {
                switch (TMath::Abs(partPdg[simuAssoc[i]]))
                {
                    case 2212:
                        protonRecEta->Fill(recoMom.PseudoRapidity());
                        protonRecP->Fill(recoMom.Mag());
                        break;
                    case 11:
                        electronRecEta->Fill(recoMom.PseudoRapidity());
                        electronRecP->Fill(recoMom.Mag());
                        break;
                    case 211:
                        pionRecEta->Fill(recoMom.PseudoRapidity());
                        pionRecP->Fill(recoMom.Mag());
                        break;
                    case 321:
                        kaonRecEta->Fill(recoMom.PseudoRapidity());
                        kaonRecP->Fill(recoMom.Mag());
                        break;
                    case 13:
                        muonRecEta->Fill(recoMom.PseudoRapidity());
                        muonRecP->Fill(recoMom.Mag());
                        break;
                    default:
                        break;
                }
            }

            initialPID->Fill(TMath::Abs(partPdg[simuAssoc[i]]));

            switch (TMath::Abs(partPdg[simuAssoc[i]]))
            {
                case 13:
                    muons[0] += 1.;
                    break;
                case 211:
                    pions[0] += 1.;
                    break;
                case 321:
                    kaons[0] += 1.;
                    break;
                case 2212:
                    protons[0] += 1.;
                    break;
                case 11:
                    electrons[0] += 1.;
                    break;
                default:
                    break;
            }

            if (trackPDG[i] == 0 || trackPDG[i] == 13)
            {
                trackPIDCutPID->Fill(TMath::Abs(partPdg[simuAssoc[i]]));
                totalParticles[1] += 1.;
                switch (TMath::Abs(partPdg[simuAssoc[i]]))
                {
                    case 13:
                        muons[1] += 1.;
                        break;
                    case 211:
                        pions[1] += 1.;
                        break;
                    case 321:
                        kaons[1] += 1.;
                        break;
                    case 2212:
                        protons[1] += 1.;
                        break;
                    case 11:
                        electrons[1] += 1.;
                        break;
                    default:
                        break;
                }
                

                if (recoMom.Mag() > 1 && std::sqrt(TMath::Abs(trackEng[i]*trackEng[i] - recoMom.Mag2())) < 0.2)
                {
                    energyCutPID->Fill(TMath::Abs(partPdg[simuAssoc[i]]));
                    totalParticles[2] += 1.;
                    switch (TMath::Abs(partPdg[simuAssoc[i]]))
                    {
                        case 13:
                            muons[2] += 1.;
                            break;
                        case 211:
                            pions[2] += 1.;
                            break;
                        case 321:
                            kaons[2] += 1.;
                            break;
                        case 2212:
                            protons[2] += 1.;
                            break;
                        case 11:
                            electrons[2] += 1.;
                            break;
                        default:
                            break;
                    }
                    bool ECalEpCut; // = (ECalEnergy/recoMom.Mag() < 0.4 && ECalEnergy/recoMom.Mag() >= 0.);
                    bool HCalEpCut;

                    if (recoMom.Mag() > 10.) // Different cuts for low-momentum tracks
                    {
                        ECalEpCut = (ECalEnergy/recoMom.Mag() <= 0.1 && ECalEnergy/recoMom.Mag() >= 0.);
                    }
                    else if (recoMom.Mag() <= 10. && recoMom.Mag() > 0.)
                    {
                        double cutValueEcal = -0.03*recoMom.Mag() + 0.4;
                        ECalEpCut = (ECalEnergy/recoMom.Mag() >= 0. && ECalEnergy/recoMom.Mag() <= cutValueEcal);
                    }
                    else
                    {
                        ECalEpCut = false;
                    }

                    if (recoMom.Mag() > 5.) // Different cuts for low-momentum tracks
                    {
                        HCalEpCut = (HCalEnergy/recoMom.Mag() <= 0.5 && HCalEnergy/recoMom.Mag() >= 0.05);
                    }
                    else if (recoMom.Mag() <= 5. && recoMom.Mag() > 0.)
                    {
                        double cutValueHcal = -0.09*recoMom.Mag() + 0.5;
                        HCalEpCut = (HCalEnergy/recoMom.Mag() >= cutValueHcal);
                    }
                    else
                    {
                        HCalEpCut = false;
                    }
                
                    if (ECalEpCut && HCalEpCut)
                    {
                        epCutPID->Fill(TMath::Abs(partPdg[simuAssoc[i]]));
                        totalParticles[3] += 1.;
                        switch (TMath::Abs(partPdg[simuAssoc[i]]))
                        {
                            case 13:
                                muons[3] += 1.;
                                break;
                            case 211:
                                pions[3] += 1.;
                                break;
                            case 321:
                                kaons[3] += 1.;
                                break;
                            case 2212:
                                protons[3] += 1.;
                                break;
                            case 11:
                                electrons[3] += 1.;
                                break;
                            default:
                                break;
                        }

                        if ((EcalHits < 4) && (HcalHits < 7))
                        {
                            clusterSizePID->Fill(TMath::Abs(partPdg[simuAssoc[i]]));
                            totalParticles[4] += 1.;
                            totalParticleRecCounts[int((TVector3(recoMom[i],recoMom[i],recoMom[i]).Mag()-0.5))] += 1;
                            switch (TMath::Abs(partPdg[simuAssoc[i]]))
                            {
                                case 13:
                                    muons[4] += 1.;
                                    muonPostCutEta->Fill(recoMom.PseudoRapidity());
                                    muonPostCutP->Fill(recoMom.Mag());
                                    muonRecCounts[int((TVector3(recoMom[i],recoMom[i],recoMom[i]).Mag()-0.5))] += 1;
                                    break;
                                case 211:
                                    pions[4] += 1.;
                                    pionPostCutEta->Fill(recoMom.PseudoRapidity());
                                    pionPostCutP->Fill(recoMom.Mag());
                                    break;
                                case 321:
                                    kaons[4] += 1.;
                                    break;
                                case 2212:
                                    protons[4] += 1.;
                                    break;
                                case 11:
                                    electrons[4] += 1.;
                                    break;
                                default:
                                    break;
                            }
                        }
                    }
                }

            }
            
        }

    }

    double sumTrueParticles = 0.;

    for (int i = 0; i < 5; i++)
    {
        sumTrueParticles += trueParticles[i];
    }

    std::cout << "" << std::endl;
    std::cout << "Event Processing Complete" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Particles Generated in Central Detector: " << sumTrueParticles << " | Initial Detected Particles: " << totalParticles[0] << " | Muons: " << muons[0] << " | Ratio: " << muons[0]/totalParticles[0] << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Efficiencies by particle type" << std::endl;
    std::cout << "Protons: " << protons[0]/trueParticles[3] << " | Electrons: " << electrons[0]/trueParticles[4] << " | Pions: " << pions[0]/trueParticles[1] << " | Kaons: " << kaons[0]/trueParticles[2] << " | Muons: " << muons[0]/trueParticles[0] << std::endl;
    std::cout << "Initital Purity of Muon Sample: " << muons[0]/totalParticles[0] <<  std::endl;
    std::cout << "Fake Rate of Muon Sample: " << (totalParticles[0]-muons[0])/totalParticles[0] <<  std::endl;
    std::cout << "" << std::endl;
    std::cout << "After Track PID Cut: " << totalParticles[1] << " | Muons: " << muons[1] << " | Ratio: " << muons[1]/totalParticles[1] << std::endl;
    std::cout << "After Energy and Mass Cut: " << totalParticles[2] << " | Muons: " << muons[2] << " | Ratio: " << muons[2]/totalParticles[2] << std::endl;
    std::cout << "After E/p Cut: " << totalParticles[3] << " | Muons: " << muons[3] << " | Ratio: " << muons[3]/totalParticles[3] << std::endl;
    std::cout << "After Cluster Size Cut: " << totalParticles[4] << " | Muons: " << muons[4] << " | Ratio: " << muons[4]/totalParticles[4] << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Efficiencies by particle type after cuts" << std::endl;
    std::cout << "Protons: " << protons[4]/trueParticles[3] << " | Electrons: " << electrons[4]/trueParticles[4] << " | Pions: " << pions[4]/trueParticles[1] << " | Kaons: " << kaons[4]/trueParticles[2] << " | Muons: " << muons[4]/trueParticles[0] << std::endl;
    std::cout << "Purity of Muon Sample after cuts: " << muons[4]/totalParticles[4] <<  std::endl;
    std::cout << "Fake Rate of Muon Sample after cuts: " << (totalParticles[4]-muons[4])/totalParticles[4] <<  std::endl;
    std::cout << "" << std::endl;

    const int n = 1000;
    double x[n+2], yE[n+2], yH[n+2];

    // line points
    for (int i = 0; i < n; ++i) {
        x[i] = 25.0 * i / (n - 1);
    }
    for (int i = 0; i < 200; ++i) {
        yH[i] = -0.09 * x[i] + 0.5;
    }
    for (int i = 200; i < n; ++i) {
        yH[i] = 0.05;
    }
    for (int i = 0; i < 400; ++i) {
        yE[i] = -0.03*x[i] + 0.4;
    }
    for (int i = 400; i < n; ++i) {
        yE[i] = 0.1;
    }

    // close the polygon down to the bottom of the histogram
    x[n]   = 25;  
    yE[n] = 2.0;
    yH[n]   = muonEpTrueEHcal->GetYaxis()->GetXmin();
    x[n+1] = 0;  
    yE[n+1] = 2.0;
    yH[n+1] = muonEpTrueEHcal->GetYaxis()->GetXmin();

    TGraph *EpEcalCuts = new TGraph(n+2, x, yE);
    EpEcalCuts->SetName("EpEcalCuts");
    EpEcalCuts->SetFillColorAlpha(kRed, 0.25);  // transparent red
    EpEcalCuts->SetLineColor(kRed);
    
    TGraph *EpHcalCuts = new TGraph(n+2, x, yH);
    EpHcalCuts->SetName("EpHcalCuts");
    EpHcalCuts->SetFillColorAlpha(kRed, 0.25);  // transparent red
    EpHcalCuts->SetLineColor(kRed);

    Int_t colours[] = {3, 2, 4, 5, 6};

    TCanvas *cEta = new TCanvas("cEta","cEta",1100,800);
    cEta->Divide(2,3);
    cEta->cd(1);
    muonRecEta->SetFillColor(colours[0]);
    muonRecEta->Draw();
    cEta->cd(2);
    pionRecEta->SetFillColor(colours[1]);
    pionRecEta->Draw();
    cEta->cd(3);
    kaonRecEta->SetFillColor(colours[2]);
    kaonRecEta->Draw();
    cEta->cd(4);
    protonRecEta->SetFillColor(colours[3]);
    protonRecEta->Draw();
    cEta->cd(5);
    electronRecEta->SetFillColor(colours[4]);
    electronRecEta->Draw();
    cEta->Update();

    TCanvas *cP = new TCanvas("cP","cP",1100,800);
    cP->Divide(2,3);
    cP->cd(1);
    muonRecP->SetFillColor(colours[0]);
    muonRecP->Draw();
    cP->cd(2);
    pionRecP->SetFillColor(colours[1]);
    pionRecP->Draw();
    cP->cd(3);
    kaonRecP->SetFillColor(colours[2]);
    kaonRecP->Draw();
    cP->cd(4);
    protonRecP->SetFillColor(colours[3]);
    protonRecP->Draw();
    cP->cd(5);
    electronRecP->SetFillColor(colours[4]);
    electronRecP->Draw();
    cP->Update();

    TCanvas *c1 = new TCanvas("c1","c1",1100,800);
    c1->Divide(2,3);
    c1->cd(1);
    muonEpTrueEcal->Draw("COLZ");
    EpEcalCuts->Draw("F SAME");
    c1->cd(2);
    pionEpTrueEcal->Draw("COLZ");
    EpEcalCuts->Draw("F SAME");
    c1->cd(3);
    kaonEpTrueEcal->Draw("COLZ");
    EpEcalCuts->Draw("F SAME");
    c1->cd(4);
    protonEpTrueEcal->Draw("COLZ");
    EpEcalCuts->Draw("F SAME");
    c1->cd(5);
    electronEpTrueEcal->Draw("COLZ");
    EpEcalCuts->Draw("F SAME");
    c1->Update();

    TCanvas *c2 = new TCanvas("c2","c2",1100,800);
    c2->Divide(2,3);
    c2->cd(1);
    muonEpRecEcal->Draw("COLZ");
    EpEcalCuts->Draw("F SAME");
    c2->cd(2);
    pionEpRecEcal->Draw("COLZ");
    EpEcalCuts->Draw("F SAME");
    c2->cd(3);
    kaonEpRecEcal->Draw("COLZ");
    EpEcalCuts->Draw("F SAME");
    c2->cd(4);
    protonEpRecEcal->Draw("COLZ");
    EpEcalCuts->Draw("F SAME");
    c2->cd(5);
    electronEpRecEcal->Draw("COLZ");
    EpEcalCuts->Draw("F SAME");
    c2->Update();

    TCanvas *c3 = new TCanvas("c3","c3",1100,800);
    c3->Divide(2,3);
    c3->cd(1);
    muonEpTrueHcal->Draw("COLZ");
    EpHcalCuts->Draw("F SAME");
    c3->cd(2);
    pionEpTrueHcal->Draw("COLZ");
    EpHcalCuts->Draw("F SAME");
    c3->cd(3);
    kaonEpTrueHcal->Draw("COLZ");
    EpHcalCuts->Draw("F SAME");
    c3->cd(4);
    protonEpTrueHcal->Draw("COLZ");
    EpHcalCuts->Draw("F SAME");
    c3->cd(5);
    electronEpTrueHcal->Draw("COLZ");
    EpHcalCuts->Draw("F SAME");
    c3->Update();

    TCanvas *c4 = new TCanvas("c4","c4",1100,800);
    c4->Divide(2,3);
    c4->cd(1);
    muonEpRecHcal->Draw("COLZ");
    EpHcalCuts->Draw("F SAME");
    c4->cd(2);
    pionEpRecHcal->Draw("COLZ");
    EpHcalCuts->Draw("F SAME");
    c4->cd(3);
    kaonEpRecHcal->Draw("COLZ");
    EpHcalCuts->Draw("F SAME");
    c4->cd(4);
    protonEpRecHcal->Draw("COLZ");
    EpHcalCuts->Draw("F SAME");
    c4->cd(5);
    electronEpRecHcal->Draw("COLZ");
    EpHcalCuts->Draw("F SAME");
    c4->Update();

    TCanvas *c5 = new TCanvas("c5","c5",1100,800);
    c5->Divide(2,3);
    c5->cd(1);
    muonEpTrueEHcal->Draw("COLZ");
    c5->cd(2);
    pionEpTrueEHcal->Draw("COLZ");
    c5->cd(3);
    kaonEpTrueEHcal->Draw("COLZ");
    c5->cd(4);
    protonEpTrueEHcal->Draw("COLZ");
    c5->cd(5);
    electronEpTrueEHcal->Draw("COLZ");
    c5->Update();

    TCanvas *c6 = new TCanvas("c6","c6",1100,800);
    c6->Divide(2,3);
    c6->cd(1);
    muonEpRecEHcal->Draw("COLZ");
    c6->cd(2);
    pionEpRecEHcal->Draw("COLZ");
    c6->cd(3);
    kaonEpRecEHcal->Draw("COLZ");
    c6->cd(4);
    protonEpRecEHcal->Draw("COLZ");
    c6->cd(5);
    electronEpRecEHcal->Draw("COLZ");
    c6->Update();

    TCanvas *c7 = new TCanvas("c7","c7",1100,800);
    c7->Divide(2,3);
    c7->cd(1);
    muonEcalClusterSize->SetFillColor(colours[0]);
    muonEcalClusterSize->Draw();
    c7->cd(2);
    pionEcalClusterSize->SetFillColor(colours[1]);
    pionEcalClusterSize->Draw();
    c7->cd(3);
    kaonEcalClusterSize->SetFillColor(colours[2]);
    kaonEcalClusterSize->Draw();
    c7->cd(4);
    protonEcalClusterSize->SetFillColor(colours[3]);
    protonEcalClusterSize->Draw();
    c7->cd(5);
    electronEcalClusterSize->SetFillColor(colours[4]);
    electronEcalClusterSize->Draw();
    c7->Update();

    TCanvas *c8 = new TCanvas("c8","c8",1100,800);
    c8->Divide(2,3);
    c8->cd(1);
    muonHcalClusterSize->SetFillColor(colours[0]);
    muonHcalClusterSize->Draw();
    c8->cd(2);
    pionHcalClusterSize->SetFillColor(colours[1]);
    pionHcalClusterSize->Draw();
    c8->cd(3);
    kaonHcalClusterSize->SetFillColor(colours[2]);
    kaonHcalClusterSize->Draw();
    c8->cd(4);
    protonHcalClusterSize->SetFillColor(colours[3]);
    protonHcalClusterSize->Draw();
    c8->cd(5);
    electronHcalClusterSize->SetFillColor(colours[4]);
    electronHcalClusterSize->Draw();
    c8->Update();

    TCanvas *c9 = new TCanvas("c9","c9",1100,800);
    c9->Divide(2,3);
    c9->cd(1);
    muonEHcalClusterSize->SetFillColor(colours[0]);
    muonEHcalClusterSize->Draw();
    c9->cd(2);
    pionEHcalClusterSize->SetFillColor(colours[1]);
    pionEHcalClusterSize->Draw();
    c9->cd(3);
    kaonEHcalClusterSize->SetFillColor(colours[2]);
    kaonEHcalClusterSize->Draw();
    c9->cd(4);
    protonEHcalClusterSize->SetFillColor(colours[3]);
    protonEHcalClusterSize->Draw();
    c9->cd(5);
    electronEHcalClusterSize->SetFillColor(colours[4]);
    electronEHcalClusterSize->Draw();
    c9->Update();


    TCanvas *cPIDcuts = new TCanvas("cPIDcuts","cPIDcuts",1100,800);
    cPIDcuts->Divide(2,3);
    cPIDcuts->cd(1);
    initialPID->Draw();
    cPIDcuts->cd(2);
    trackPIDCutPID->Draw();
    cPIDcuts->cd(3);
    energyCutPID->Draw();
    cPIDcuts->cd(4);
    epCutPID->Draw();
    cPIDcuts->cd(5);
    clusterSizePID->Draw();
    cPIDcuts->Update();

    for (size_t i = 0; i < 20; i++)
    {
        if (totalParticleTrueCounts[i] == 0) totalParticleTrueCounts[i] = 1;
        if (totalParticleRecCounts[i] == 0) totalParticleRecCounts[i] = 1; 
        muonPurityTrue[i] = muonTrueCounts[i]/(totalParticleTrueCounts[i]/2.4);   
        muonPurityRec[i] = muonRecCounts[i]/totalParticleRecCounts[i];
        muonEfficiency[i] = muonRecCounts[i]/muonTrueCounts[i];
    }

    TGraph *muonEff = new TGraph(20, momentumBins, muonEfficiency);
    muonEff->SetTitle("Muon Efficiency vs Momentum;Momentum (GeV/c);Efficiency");
    TGraph *muonPurR = new TGraph(20, momentumBins, muonPurityRec);
    muonPurR->SetTitle("Muon Purity vs Momentum;Momentum (GeV/c);Purity");
    TGraph *muonPurT = new TGraph(20, momentumBins, muonPurityTrue);
    muonPurT->SetTitle("Muon Purity vs Momentum;Momentum (GeV/c);Purity");

    TCanvas *cMuonEffPur = new TCanvas("cMuonEffPur","cMuonEffPur",1100,800);
    cMuonEffPur->Divide(2,2);
    cMuonEffPur->cd(1);
    muonEff->SetLineColor(kBlue);
    muonEff->SetMarkerColor(kBlue);
    muonEff->SetMarkerStyle(20);
    muonEff->Draw("APL");
    cMuonEffPur->cd(2);
    muonPurT->SetLineColor(kRed);
    muonPurT->SetMarkerColor(kRed);
    muonPurT->SetMarkerStyle(20);
    muonPurT->Draw("APL");
    cMuonEffPur->cd(3);
    muonPurR->SetLineColor(kGreen+2);
    muonPurR->SetMarkerColor(kGreen+2);
    muonPurR->SetMarkerStyle(20);
    muonPurR->Draw("APL");
    cMuonEffPur->Update();


    for (int i = 0; i < 5; i++)
    {
        if (muons[i] == 0) muons[i] = 1;
        if (pions[i] == 0) pions[i] = 1;
        if (kaons[i] == 0) kaons[i] = 1;
        if (protons[i] == 0) protons[i] = 1;
        if (electrons[i] == 0) electrons[i] = 1;
    }

    TCanvas *pieCharts = new TCanvas("pieCharts","pieCharts",1100,800);
    pieCharts->Divide(2,3);
    pieCharts->cd(1);
    TPie *initialPIDPie = new TPie("initialPIDPie", "Initial Particle ID Distribution", 5);
    initialPIDPie->SetEntryVal(0, muons[0]);
    initialPIDPie->SetEntryVal(1, pions[0]); 
    initialPIDPie->SetEntryVal(2, kaons[0]); 
    initialPIDPie->SetEntryVal(3, protons[0]); 
    initialPIDPie->SetEntryVal(4, electrons[0]); 
    initialPIDPie->SetEntryLabel(0, "Muons");
    initialPIDPie->SetEntryLabel(1, "Pions");
    initialPIDPie->SetEntryLabel(2, "Kaons");
    initialPIDPie->SetEntryLabel(3, "Protons");
    initialPIDPie->SetEntryLabel(4, "Electrons");
    for (int i = 0; i < 5; i++)
    {
        initialPIDPie->SetEntryFillColor(i, colours[i]);
    }
    initialPIDPie->SetTitle("Initial Particle ID Distribution");
    initialPIDPie->Draw();
    pieCharts->cd(2);
    TPie *trackPIDCutPie = new TPie("trackPIDCutPie", "Particle ID Distribution After Track PID Cut", 5);
    trackPIDCutPie->SetEntryVal(0, muons[1]);  
    trackPIDCutPie->SetEntryVal(1, pions[1]); 
    trackPIDCutPie->SetEntryVal(2, kaons[1]); 
    trackPIDCutPie->SetEntryVal(3, protons[1]); 
    trackPIDCutPie->SetEntryVal(4, electrons[1]);
    trackPIDCutPie->SetEntryLabel(0,"Muons");
    trackPIDCutPie->SetEntryLabel(1,"Pions");
    trackPIDCutPie->SetEntryLabel(2,"Kaons");
    trackPIDCutPie->SetEntryLabel(3,"Protons");
    trackPIDCutPie->SetEntryLabel(4,"Electrons");
    for (int i = 0; i < 5; i++)
    {
        trackPIDCutPie->SetEntryFillColor(i, colours[i]);
    }
    trackPIDCutPie->SetTitle("After Track PID Cut");
    trackPIDCutPie->Draw();
    pieCharts->cd(3);
    TPie *energyCutPie = new TPie("energyCutPie", "Particle ID Distribution After Energy Cut", 5);
    energyCutPie->SetEntryVal(0, muons[2]);  
    energyCutPie->SetEntryVal(1, pions[2]); 
    energyCutPie->SetEntryVal(2, kaons[2]); 
    energyCutPie->SetEntryVal(3, protons[2]); 
    energyCutPie->SetEntryVal(4, electrons[2]);
    energyCutPie->SetEntryLabel(0,"Muons");
    energyCutPie->SetEntryLabel(1,"Pions");
    energyCutPie->SetEntryLabel(2,"Kaons");
    energyCutPie->SetEntryLabel(3,"Protons");
    energyCutPie->SetEntryLabel(4,"Electrons");
    for (int i = 0; i < 5; i++)
    {
        energyCutPie->SetEntryFillColor(i, colours[i]);
    }
    energyCutPie->SetTitle("After Energy Cut");
    energyCutPie->Draw();
    pieCharts->cd(4);
    TPie *epCutPie = new TPie("epCutPie", "Particle ID Distribution After E/p Cut", 5);
    epCutPie->SetEntryVal(0, muons[3]);  
    epCutPie->SetEntryVal(1, pions[3]); 
    epCutPie->SetEntryVal(2, kaons[3]); 
    epCutPie->SetEntryVal(3, protons[3]); 
    epCutPie->SetEntryVal(4, electrons[3]);
    epCutPie->SetEntryLabel(0,"Muons");
    epCutPie->SetEntryLabel(1,"Pions");
    epCutPie->SetEntryLabel(2,"Kaons");
    epCutPie->SetEntryLabel(3,"Protons");
    epCutPie->SetEntryLabel(4,"Electrons");
    for (int i = 0; i < 5; i++)
    {
        epCutPie->SetEntryFillColor(i, colours[i]);
    }
    epCutPie->SetTitle("After E/p Cut");
    epCutPie->Draw();
    pieCharts->cd(5);
    TPie *clusterSizePie = new TPie("clusterSizePie", "Particle ID Distribution After Cluster Size Cut", 5);
    clusterSizePie->SetEntryVal(0, muons[4]);  
    clusterSizePie->SetEntryVal(1, pions[4]); 
    clusterSizePie->SetEntryVal(2, kaons[4]); 
    clusterSizePie->SetEntryVal(3, protons[4]); 
    clusterSizePie->SetEntryVal(4, electrons[4]);
    clusterSizePie->SetEntryLabel(0,"Muons");
    clusterSizePie->SetEntryLabel(1,"Pions");
    clusterSizePie->SetEntryLabel(2,"Kaons");
    clusterSizePie->SetEntryLabel(3,"Protons");
    clusterSizePie->SetEntryLabel(4,"Electrons");
    for (int i = 0; i < 5; i++)
    {
        clusterSizePie->SetEntryFillColor(i, colours[i]);
    }
    clusterSizePie->SetTitle("After Cluster Size Cut");
    clusterSizePie->Draw();
    pieCharts->Update();

    TCanvas *cPostCuts = new TCanvas("cPostCuts","cPostCuts",1100,800);
    cPostCuts->Divide(2,3);
    cPostCuts->cd(1);
    muonPostCutEta->SetFillColor(colours[0]);
    muonPostCutEta->Draw();
    cPostCuts->cd(2);
    muonPostCutP->SetFillColor(colours[0]);
    muonPostCutP->Draw();
    cPostCuts->cd(3);
    pionPostCutEta->SetFillColor(colours[1]);
    pionPostCutEta->Draw();
    cPostCuts->cd(4);
    pionPostCutP->SetFillColor(colours[1]);
    pionPostCutP->Draw();
    cPostCuts->Update();

    ofile->cd();
    ofile->mkdir("eta");
    ofile->cd("eta");
    protonRecEta->Write();
    electronRecEta->Write();
    pionRecEta->Write();
    kaonRecEta->Write();
    muonRecEta->Write();
    ofile->cd("..");
    ofile->mkdir("momentum");
    ofile->cd("momentum");
    protonRecP->Write();
    electronRecP->Write();
    pionRecP->Write();
    kaonRecP->Write();
    muonRecP->Write();
    ofile->cd("..");
    ofile->mkdir("trueEcalEp");
    ofile->cd("trueEcalEp");
    EpEcalCuts->Write();
    protonEpTrueEcal->Write();
    electronEpTrueEcal->Write();
    pionEpTrueEcal->Write();
    kaonEpTrueEcal->Write();
    muonEpTrueEcal->Write();
    ofile->cd("..");
    ofile->mkdir("recoEcalEp");
    ofile->cd("recoEcalEp");
    EpEcalCuts->Write();
    protonEpRecEcal->Write();
    electronEpRecEcal->Write();
    pionEpRecEcal->Write();
    kaonEpRecEcal->Write();
    muonEpRecEcal->Write();
    ofile->cd("..");
    ofile->mkdir("trueHcalEp");
    ofile->cd("trueHcalEp");
    EpHcalCuts->Write();
    protonEpTrueHcal->Write();
    electronEpTrueHcal->Write();
    pionEpTrueHcal->Write();
    kaonEpTrueHcal->Write();
    muonEpTrueHcal->Write();
    ofile->cd("..");
    ofile->mkdir("recoHcalEp");
    ofile->cd("recoHcalEp");
    EpHcalCuts->Write();
    protonEpRecHcal->Write();
    electronEpRecHcal->Write();
    pionEpRecHcal->Write();
    kaonEpRecHcal->Write();
    muonEpRecHcal->Write();
    ofile->cd("..");
    ofile->mkdir("trueEHcal");
    ofile->cd("trueEHcal");
    protonEpTrueEHcal->Write();
    electronEpTrueEHcal->Write();
    pionEpTrueEHcal->Write();
    kaonEpTrueEHcal->Write();
    muonEpTrueEHcal->Write();
    ofile->cd("..");
    ofile->mkdir("recoEHcal");
    ofile->cd("recoEHcal");
    protonEpRecEHcal->Write();
    electronEpRecEHcal->Write();
    pionEpRecEHcal->Write();
    kaonEpRecEHcal->Write();
    muonEpRecEHcal->Write();
    ofile->cd("..");
    ofile->mkdir("CalClusterSize");
    ofile->cd("CalClusterSize");
    protonEcalClusterSize->Write();
    electronEcalClusterSize->Write();
    pionEcalClusterSize->Write();
    kaonEcalClusterSize->Write();
    muonEcalClusterSize->Write();
    protonHcalClusterSize->Write();
    electronHcalClusterSize->Write();
    pionHcalClusterSize->Write();
    kaonHcalClusterSize->Write();
    muonHcalClusterSize->Write();
    protonEHcalClusterSize->Write();
    electronEHcalClusterSize->Write();
    pionEHcalClusterSize->Write();
    kaonEHcalClusterSize->Write();
    muonEHcalClusterSize->Write();
    ofile->cd("..");
    ofile->mkdir("cuts");
    ofile->cd("cuts");
    initialPID->Write();
    initialPIDPie->Write();
    trackPIDCutPID->Write();
    trackPIDCutPie->Write();
    energyCutPID->Write();
    energyCutPie->Write();
    epCutPID->Write();
    epCutPie->Write();
    clusterSizePID->Write();
    clusterSizePie->Write();
    ofile->cd("..");
    ofile->mkdir("postCuts");
    ofile->cd("postCuts");
    muonPostCutEta->Write();
    muonPostCutP->Write();
    pionPostCutEta->Write();
    pionPostCutP->Write();
    muonEff->Write("muonEfficiency");
    muonPurT->Write("muonPurityTrue");
    muonPurR->Write("muonPurityRec");
    ofile->cd("..");


}