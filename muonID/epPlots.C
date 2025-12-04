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
#include <TMultiGraph.h>
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

#include "../ePICStyle.C"

void epPlots()
{
    //gROOT->SetBatch(kTRUE);
    gROOT->ProcessLine("SetePICStyle()");
    gStyle->SetOptStat(0);

    //TString infile="../eicReconOutput/EICreconOut_JPsiMuMu_10ifb_10x130ep_Pruned.root";
    TString infile="reconOut/mu-pi_20GeV_reconOut.root";
    //TString infile="../dis_background/DIS_Q2_1_10_10x130ep_Pruned.root";

    //std::string outfilename = "outputs/DIS_Q2_1_10_10x130ep_epPlots.root";
    std::string outfilename = "outputs/mu-pi_epPlots.root";
    //std::string outfilename = "outputs/JPsiMuMu_10ifb_10x130ep_epPlots.root";

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

    TTreeReaderArray<unsigned int> simuAssocLFHcal(tree_reader, "LFHCALClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocLFHcal(tree_reader, "LFHCALClusterAssociations.recID");    
    TTreeReaderArray<float> LFHcalEng(tree_reader, "LFHCALClusters.energy");

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
    TH1D *allParticlesTrueP = new TH1D("allParticlesTrueP","Momentum of all Generated Particles; p [GeV/c]; Counts",50,0.,25.);
    TH1D *allParticlesTrueP_PC = new TH1D("allParticlesTrueP_PC","Momentum of all Reconstructed Particles post cuts; p [GeV/c]; Counts",50,0.,25.);

    TH2D *allParticlesTrueEtaP =  new TH2D("allParticlesTrueEtaP", "Momentum of all Generated Particles against Eta; p [GeV/c]; eta",50,0.,25.,100,-5.,5);
    TH2D *allParticlesTrueEtaP_PC =  new TH2D("allParticlesTrueEtaP_PC", "Momentum of all Reconstructed Particles against Eta post cuts; p [GeV/c]; eta",50,0.,25.,100,-5.,5);

    TH1D *protonRecP = new TH1D("protonP","Reconstructed Proton Momentum; p [GeV/c]; Counts",50,0.,25.);
    TH1D *electronRecP = new TH1D("electronP","Reconstructed Electron Momentum; p [GeV/c]; Counts",50,0.,25.);
    TH1D *pionRecP = new TH1D("pionP","Reconstructed Pion Momentum; p [GeV/c]; Counts",50,0.,25.);
    TH1D *kaonRecP = new TH1D("kaonP","Reconstructed Kaon Momentum; p [GeV/c]; Counts",50,0.,25.);

    TH1D *muonTrueP = new TH1D("muonTrueP","Momentum of Generated Muons;p [GeV/c]",50,0.,25.);
    TH2D *muonTruePEta = new TH2D("muonTruePEta","Momentum of Generated Muons against eta; p [GeV/c]; eta",50,0.,25.,100,-5.,5);

    TH2D *muonRecPEta = new TH2D("muonRecPEta","Reconstructed of Generated Muons against eta; p [GeV/c]; eta",50,0.,25.,100,-5.,5);
    TH1D *muonRecP = new TH1D("muonP","Reconstructed Muon Momentum; p [GeV/c]; Counts",50,0.,25.);
    
    TH1D *muonMomEff = new TH1D("muonMomEff","Efficency;Momentum (GeV/c)",50,0.,25.);
    TH1D *muonMomPur = new TH1D("muonMomPur","Purity;Momentum (GeV/c)",50,0.,25.);

    TH2D *muonMomEtaEff = new TH2D("muonMomEtaEff","Efficency as a function of eta and momentum; p [GeV/c]; eta",50,0.,25.,100,-5.,5);
    TH2D *muonMomEtaPur = new TH2D("muonMomEtaPur","Purity as a function of eta and momentum; p [GeV/c]; eta",50,0.,25.,100,-5.,5);

    TH1D *muonRecMinusTrueP = new TH1D("muonRecMinusTrueP","Muon Momentum Difference; p [GeV/c]; Counts",100,-25.,25.);

    // Ecal ep Plots by particle
    TH2D *protonEpTrueEcal = new TH2D("protonEpTrueEcal","Ecal Energy vs Track Momentum for Protons; p;  E/p",50,0.,25.,200,0.,2.0);
    TH2D *electronEpTrueEcal = new TH2D("electronEpTrueEcal","Ecal Energy vs Track Momentum for Electrons; p;  E/p",50,0.,25.,100,0.,2.0);
    TH2D *pionEpTrueEcal = new TH2D("pionEpTrueEcal","Ecal Energy vs Track Momentum for Pions; p;  E/p",50,0.,25.,100,0.,2.0);
    TH2D *kaonEpTrueEcal = new TH2D("kaonEpTrueEcal","Ecal Energy vs Track Momentum for Kaons; p;  E/p",50,0.,25.,100,0.,2.0);
    TH2D *muonEpTrueEcal = new TH2D("muonEpTrueEcal","Ecal Energy vs Track Momentum for Muons; p;  E/p",50,0.,25.,100,0.,2.0);

    TH2D *protonEpRecEcal = new TH2D("protonEpRecEcal","Ecal Energy vs Track Momentum for Protons; p;  E/p",50,0.,25.,200,0.,2.0);
    TH2D *electronEpRecEcal = new TH2D("electronEpRecEcal","Ecal Energy vs Track Momentum for Electrons; p;  E/p",50,0.,25.,100,0.,2.0);
    TH2D *pionEpRecEcal = new TH2D("pionEpRecEcal","Ecal Energy vs Track Momentum for Pions; p;  E/p",50,0.,25.,100,0.,2.0);
    TH2D *kaonEpRecEcal = new TH2D("kaonEpRecEcal","Ecal Energy vs Track Momentum for Kaons; p;  E/p",50,0.,25.,100,0.,2.0);
    TH2D *muonEpRecEcal = new TH2D("muonEpRecEcal","Ecal Energy vs Track Momentum for Muons; p;  E/p",50,0.,25.,100,0.,2.0);

    TH2D *mupiEpTrueEcalvsEta = new TH2D("mupiEpTrueEcalvsEta","Ecal Energy vs Track Eta for all particles; eta;  E/p",100,-5.,5.,100,0.,2.0);
    TH2D *mupiEpRecEcalvsEta = new TH2D("mupiEpRecEcalvsEta","Ecal Energy vs Track Eta for all particles; eta;  E/p",100,-5.,5.,100,0.,2.0);

    // Hcal ep Plots by particle
    TH2D *protonEpTrueHcal = new TH2D("protonEpTrueHcal","Hcal Energy vs Track Momentum for Protons; p;  E/p",50,0.,25.,200,0.,2.0);
    TH2D *electronEpTrueHcal = new TH2D("electronEpTrueHcal","Hcal Energy vs Track Momentum for Electrons; p;  E/p",50,0.,25.,100,0.,2.0);
    TH2D *pionEpTrueHcal = new TH2D("pionEpTrueHcal","Hcal Energy vs Track Momentum for Pions; p;  E/p",50,0.,25.,100,0.,2.0);
    TH2D *kaonEpTrueHcal = new TH2D("kaonEpTrueHcal","Hcal Energy vs Track Momentum for Kaons; p;  E/p",50,0.,25.,100,0.,2.0);
    TH2D *muonEpTrueHcal = new TH2D("muonEpTrueHcal","Hcal Energy vs Track Momentum for Muons; p;  E/p",50,0.,25.,100,0.,2.0);

    TH2D *protonEpRecHcal = new TH2D("protonEpRecHcal","Hcal Energy vs Track Momentum for Protons; p;  E/p",50,0.,25.,200,0.,2.0);
    TH2D *electronEpRecHcal = new TH2D("electronEpRecHcal","Hcal Energy vs Track Momentum for Electrons; p;  E/p",50,0.,25.,100,0.,2.0);
    TH2D *pionEpRecHcal = new TH2D("pionEpRecHcal","Hcal Energy vs Track Momentum for Pions; p;  E/p",50,0.,25.,100,0.,2.0);
    TH2D *kaonEpRecHcal = new TH2D("kaonEpRecHcal","Hcal Energy vs Track Momentum for Kaons; p;  E/p",50,0.,25.,100,0.,2.0);
    TH2D *muonEpRecHcal = new TH2D("muonEpRecHcal","Hcal Energy vs Track Momentum for Muons; p;  E/p",50,0.,25.,100,0.,2.0);

    TH2D *mupiEpTrueHcalvsEta = new TH2D("mupiEpTrueHcalvsEta","Hcal Energy vs Track Eta for all particles; eta;  E/p",100,-5.,5.,100,0.,2.0);
    TH2D *mupiEpRecHcalvsEta = new TH2D("mupiEpRecHcalvsEta","Hcal Energy vs Track Eta for all particles; eta;  E/p",100,-5.,5.,100,0.,2.0);

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
    TH2D *protonEcalClusterSize = new TH2D("protonEcalClusterSize","Ecal Barrel Cluster Size for Protons; Momentum (GeV/c); Size",50,0.,25.,100,0.,50.);
    TH2D *electronEcalClusterSize = new TH2D("electronEcalClusterSize","Ecal Barrel Cluster Size for Electrons; Momentum (GeV/c); Size",50,0.,25.,100,0.,50.);
    TH2D *pionEcalClusterSize = new TH2D("pionEcalClusterSize","Ecal Barrel Cluster Size for Pions; Momentum (GeV/c); Size",50,0.,25.,100,0.,50.);
    TH2D *kaonEcalClusterSize = new TH2D("kaonEcalClusterSize","Ecal Barrel Cluster Size for Kaons; Momentum (GeV/c); Size",50,0.,25.,100,0.,50.);
    TH2D *muonEcalClusterSize = new TH2D("muonEcalClusterSize","Ecal Barrel Cluster Size for Muons; Momentum (GeV/c); Size",50,0.,25.,100,0.,50.);

    TH2D *protonHcalClusterSize = new TH2D("protonHcalClusterSize","Hcal Barrel Cluster Size for Protons; Momentum (GeV/c); Size",50,0.,25.,100,0.,50.);
    TH2D *electronHcalClusterSize = new TH2D("electronHcalClusterSize","Hcal Barrel Cluster Size for Electrons; Momentum (GeV/c); Size",50,0.,25.,100,0.,50.);
    TH2D *pionHcalClusterSize = new TH2D("pionHcalClusterSize","Hcal Barrel Cluster Size for Pions; Momentum (GeV/c); Size",50,0.,25.,100,0.,50.);
    TH2D *kaonHcalClusterSize = new TH2D("kaonHcalClusterSize","Hcal Barrel Cluster Size for Kaons; Momentum (GeV/c); Size",50,0.,25.,100,0.,50.);
    TH2D *muonHcalClusterSize = new TH2D("muonHcalClusterSize","Hcal Barrel Cluster Size for Muons; Momentum (GeV/c); Size",50,0.,25.,100,0.,50.);

    TH2D *protonEHcalClusterSize = new TH2D("protonEHcalClusterSize","Ecal + Hcal Barrel Cluster Size for Protons; Momentum (GeV/c); Size",50,0.,25.,200,0.,100.);
    TH2D *electronEHcalClusterSize = new TH2D("electronEHcalClusterSize","Ecal + Hcal Barrel Cluster Size for Electrons; Momentum (GeV/c); Size",50,0.,25.,200,0.,100.);
    TH2D *pionEHcalClusterSize = new TH2D("pionEHcalClusterSize","Ecal + Hcal Barrel Cluster Size for Pions; Momentum (GeV/c); Size",50,0.,25.,200,0.,100.);
    TH2D *kaonEHcalClusterSize = new TH2D("kaonEHcalClusterSize","Ecal + Hcal Barrel Cluster Size for Kaons; Momentum (GeV/c); Size",50,0.,25.,200,0.,100.);
    TH2D *muonEHcalClusterSize = new TH2D("muonEHcalClusterSize","Ecal + Hcal Barrel Cluster Size for Muons; Momentum (GeV/c); Size",50,0.,25.,200,0.,100.);

    // Post cut plot

    TH1D *initialPID = new TH1D("initialPID","Initial Particle ID Distribution; PDG; Counts",2500,0,2500);
    TH1D *trackPIDCutPID = new TH1D("trackPIDCutPID","Particle ID Distribution after Track PID Cut; PDG; Counts",2500,0,2500);
    TH1D *energyCutPID = new TH1D("energyCutPID","Particle ID Distribution after Energy and Mass Cut; PDG; Counts",2500,0,2500);
    TH1D *epCutPID = new TH1D("epCutPID","Particle ID Distribution after E/p Cut; PDG; Counts",2500,0,2500);
    TH1D *clusterSizePID = new TH1D("clusterSizePID","Particle ID Distribution after Cluster Size Cut; PDG; Counts",2500,0,2500);

    TH1D *pionPostCutEta = new TH1D("pionPCEta","Pion Pseudo-rapidity post cuts; eta; Counts",100,-5.,5.);
    TH1D *muonPostCutEta = new TH1D("muonPCEta","Muon Pseudo-rapidity post cuts; eta; Counts",100,-5.,5.);

    TH1D *pionPostCutP = new TH1D("pionPCP","Pion Momentum post cuts; p [GeV/c]; Counts",50,0.,25.);
    TH1D *muonPostCutP = new TH1D("muonPCP","Muon Momentum post cuts; p [GeV/c]; Counts",50,0.,25.);

    // E/p cut functions

    TF1 *eCalCutFunction = new TF1("eCalCutFunction","[0]/(x + [1])",0.,20.,0.,20.);
    eCalCutFunction->SetParameter(0,0.29);
    eCalCutFunction->SetParameter(1,0.23);

    TF1 *eCalCutWidthFunction = new TF1("eCalCutWidthFunction","3*([0]/x)",0.,20.);
    eCalCutWidthFunction->SetParameter(0,0.09);

    TF1 *hCalCutFunction = new TF1("hCalCutFunction","[0]/(x + [1])",0.,20.,0.,20.);
    hCalCutFunction->SetParameter(0,1.03);
    hCalCutFunction->SetParameter(1,0.44);

    TF1 *hCalCutWidthFunction = new TF1("hCalCutWidthFunction","3*([0]/x)",0.,20.);
    hCalCutWidthFunction->SetParameter(0,0.21);

    // Define 3 and 4-momentum vectors for particles

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

    bool piEcalSlices, muEcalSlices, piHcalSlices, muHcalSlices = false;

    while(tree_reader.Next()) // Loop over events
    {
      eventID++;

      if (eventID % 10 == 0)
      {
        fprintf (stderr, "%4.2f Percent\r ", eventID*100.0/mychain->GetEntries());
        fflush (stderr);
      }

      // Reset vectors
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

            if (partEta < -4 || partEta > 4) continue;
            if (TVector3(partMomX[i],partMomY[i],partMomZ[i]).Mag() < 0.05 || TVector3(partMomX[i],partMomY[i],partMomZ[i]).Mag() > 20) continue;
          
            allParticlesTrueP->Fill(TVector3(partMomX[i],partMomY[i],partMomZ[i]).Mag());
            allParticlesTrueEtaP->Fill(TVector3(partMomX[i],partMomY[i],partMomZ[i]).Mag(), TVector3(partMomX[i],partMomY[i],partMomZ[i]).Eta());

            switch (TMath::Abs(partPdg[i]))
            {
                case 13:
                    trueParticles[0] += 1.;
                    muonTrueP->Fill(TVector3(partMomX[i],partMomY[i],partMomZ[i]).Mag());
                    muonTruePEta->Fill(TVector3(partMomX[i],partMomY[i],partMomZ[i]).Mag(), TVector3(partMomX[i],partMomY[i],partMomZ[i]).Eta());
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
            if (recoMom.Mag() < 0.05 || recoMom.Mag() > 20.0) continue;

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

            if (LFHcalEng.GetSize() == simuAssocLFHcal.GetSize())
            {
                for (int jL = 0; jL < simuAssocLFHcal.GetSize(); jL++) // Look for associations in the LFHcal
                {
                    if (simuAssocLFHcal[jL] == simuAssoc[i])
                    {
                        HCalEnergy += LFHcalEng[jL];
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
                        protonEcalClusterSize->Fill(recoMom.Mag(),EcalHits);
                        protonHcalClusterSize->Fill(recoMom.Mag(),HcalHits);
                        protonEHcalClusterSize->Fill(recoMom.Mag(),EcalHits+HcalHits);
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
                        electronEcalClusterSize->Fill(recoMom.Mag(),EcalHits);
                        electronHcalClusterSize->Fill(recoMom.Mag(),HcalHits);
                        electronEHcalClusterSize->Fill(recoMom.Mag(),EcalHits+HcalHits);
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
                        pionEcalClusterSize->Fill(recoMom.Mag(),EcalHits);
                        pionHcalClusterSize->Fill(recoMom.Mag(),HcalHits);
                        pionEHcalClusterSize->Fill(recoMom.Mag(),EcalHits+HcalHits);
                        mupiEpRecEcalvsEta->Fill(recoMom.PseudoRapidity(), ECalEnergy/recoMom.Mag());
                        mupiEpTrueEcalvsEta->Fill(trueMom.PseudoRapidity(), ECalEnergy/trueMom.Mag());
                        mupiEpRecHcalvsEta->Fill(recoMom.PseudoRapidity(), HCalEnergy/recoMom.Mag());
                        mupiEpTrueHcalvsEta->Fill(trueMom.PseudoRapidity(), HCalEnergy/trueMom.Mag());
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
                        kaonEcalClusterSize->Fill(recoMom.Mag(),EcalHits);
                        kaonHcalClusterSize->Fill(recoMom.Mag(),HcalHits);
                        kaonEHcalClusterSize->Fill(recoMom.Mag(),EcalHits+HcalHits);
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
                        muonEcalClusterSize->Fill(recoMom.Mag(),EcalHits);
                        muonHcalClusterSize->Fill(recoMom.Mag(),HcalHits);
                        muonEHcalClusterSize->Fill(recoMom.Mag(),EcalHits+HcalHits);
                        mupiEpRecEcalvsEta->Fill(recoMom.PseudoRapidity(), ECalEnergy/recoMom.Mag());
                        mupiEpTrueEcalvsEta->Fill(trueMom.PseudoRapidity(), ECalEnergy/trueMom.Mag());
                        mupiEpRecHcalvsEta->Fill(recoMom.PseudoRapidity(), HCalEnergy/recoMom.Mag());
                        mupiEpTrueHcalvsEta->Fill(trueMom.PseudoRapidity(), HCalEnergy/trueMom.Mag());
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
                        protonEcalClusterSize->Fill(recoMom.Mag(),EcalHits);
                        break;
                    case 11:
                        electronRecEta->Fill(recoMom.PseudoRapidity());
                        electronRecP->Fill(recoMom.Mag());
                        electronEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        electronEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        electronEcalClusterSize->Fill(recoMom.Mag(),EcalHits);
                        break;
                    case 211:
                        pionRecEta->Fill(recoMom.PseudoRapidity());
                        pionRecP->Fill(recoMom.Mag());
                        pionEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        pionEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        pionEcalClusterSize->Fill(recoMom.Mag(),EcalHits);
                        mupiEpRecEcalvsEta->Fill(recoMom.PseudoRapidity(), ECalEnergy/recoMom.Mag());
                        mupiEpTrueEcalvsEta->Fill(trueMom.PseudoRapidity(), ECalEnergy/trueMom.Mag());
                        break;
                    case 321:
                        kaonRecEta->Fill(recoMom.PseudoRapidity());
                        kaonRecP->Fill(recoMom.Mag());
                        kaonEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        kaonEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        kaonEcalClusterSize->Fill(recoMom.Mag(),EcalHits);
                        break;
                    case 13:
                        muonRecEta->Fill(recoMom.PseudoRapidity());
                        muonRecP->Fill(recoMom.Mag());
                        muonEpTrueEcal->Fill(trueMom.Mag(), ECalEnergy/trueMom.Mag());
                        muonEpRecEcal->Fill(recoMom.Mag(), ECalEnergy/recoMom.Mag());
                        muonEcalClusterSize->Fill(recoMom.Mag(),EcalHits);
                        mupiEpRecEcalvsEta->Fill(recoMom.PseudoRapidity(), ECalEnergy/recoMom.Mag());
                        mupiEpTrueEcalvsEta->Fill(trueMom.PseudoRapidity(), ECalEnergy/trueMom.Mag());
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
                        protonHcalClusterSize->Fill(recoMom.Mag(),HcalHits);
                        break;
                    case 11:
                        electronRecEta->Fill(recoMom.PseudoRapidity());
                        electronRecP->Fill(recoMom.Mag());
                        electronEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        electronEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        electronHcalClusterSize->Fill(recoMom.Mag(),HcalHits);
                        break;
                    case 211:
                        pionRecEta->Fill(recoMom.PseudoRapidity());
                        pionRecP->Fill(recoMom.Mag());
                        pionEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        pionEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        pionHcalClusterSize->Fill(recoMom.Mag(),HcalHits);
                        mupiEpRecHcalvsEta->Fill(recoMom.PseudoRapidity(), HCalEnergy/recoMom.Mag());
                        mupiEpTrueHcalvsEta->Fill(trueMom.PseudoRapidity(), HCalEnergy/trueMom.Mag());
                        break;
                    case 321:
                        kaonRecEta->Fill(recoMom.PseudoRapidity());
                        kaonRecP->Fill(recoMom.Mag());
                        kaonEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        kaonEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        kaonHcalClusterSize->Fill(recoMom.Mag(),HcalHits);
                        break;
                    case 13:
                        muonRecEta->Fill(recoMom.PseudoRapidity());
                        muonRecP->Fill(recoMom.Mag());
                        muonEpTrueHcal->Fill(trueMom.Mag(), HCalEnergy/trueMom.Mag());
                        muonEpRecHcal->Fill(recoMom.Mag(), HCalEnergy/recoMom.Mag());
                        muonHcalClusterSize->Fill(recoMom.Mag(),HcalHits);        
                        mupiEpRecHcalvsEta->Fill(recoMom.PseudoRapidity(), HCalEnergy/recoMom.Mag());
                        mupiEpTrueHcalvsEta->Fill(trueMom.PseudoRapidity(), HCalEnergy/trueMom.Mag());
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

                    if (recoMom.Mag() > 0.5) // Different cuts for low-momentum tracks
                    {
                        double cutValueEcal = eCalCutFunction->Eval(recoMom.Mag());
                        double cutWidthEcal = eCalCutWidthFunction->Eval(recoMom.Mag());

                        ECalEpCut = (ECalEnergy/recoMom.Mag() >= 0.0 && ECalEnergy/recoMom.Mag() <= cutValueEcal+cutWidthEcal);
                    }
                    else
                    {
                        ECalEpCut = false;
                    }

                    if (recoMom.Mag() > 0.5) // Different cuts for low-momentum tracks
                    {
                        double cutValueHcal = hCalCutFunction->Eval(recoMom.Mag());
                        double cutWidthHcal = hCalCutWidthFunction->Eval(recoMom.Mag());

                        HCalEpCut = (HCalEnergy/recoMom.Mag() >= cutValueHcal-cutWidthHcal && HCalEnergy/recoMom.Mag() <= cutValueHcal+cutWidthHcal);
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

                        if ((EcalHits < 5) && (HcalHits < 8))
                        {
                            clusterSizePID->Fill(TMath::Abs(partPdg[simuAssoc[i]]));
                            allParticlesTrueP_PC->Fill(trueMom.Mag());
                            allParticlesTrueEtaP_PC->Fill(trueMom.Mag(), trueMom.Eta());
                            totalParticles[4] += 1.;
                            switch (TMath::Abs(partPdg[simuAssoc[i]]))
                            {
                                case 13:
                                    muons[4] += 1.;
                                    muonPostCutEta->Fill(recoMom.PseudoRapidity());
                                    muonPostCutP->Fill(recoMom.Mag());
                                    muonRecMinusTrueP->Fill(recoMom.Mag()-trueMom.Mag());
                                    muonMomEff->Fill(trueMom.Mag());
                                    muonMomEtaEff->Fill(trueMom.Mag(), trueMom.Eta());
                                    muonMomPur->Fill(trueMom.Mag());
                                    muonMomEtaPur->Fill(trueMom.Mag(), trueMom.Eta());
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
    double x[n], yE[n], yH[n];
    double yEup[n], yHup[n], yEdn[n], yHdn[n];

    // line points
    for (int i = 0; i < n; ++i) {
        x[i] = (19.5 * i / (n - 1)) + 0.5; // from 0.5 to 20 GeV
    }
    for (int i = 0; i < n; ++i) {
        yE[i] = eCalCutFunction->Eval(x[i]);
        yEup[i] = eCalCutFunction->Eval(x[i]) + eCalCutWidthFunction->Eval(x[i]);
        yEdn[i] = 0.0;
    }
    for (int i = 0; i < n; ++i) {
        yH[i] = hCalCutFunction->Eval(x[i]);
        yHup[i] = hCalCutFunction->Eval(x[i]) + hCalCutWidthFunction->Eval(x[i]);
        yHdn[i] = hCalCutFunction->Eval(x[i]) - hCalCutWidthFunction->Eval(x[i]);
    }

    TGraph *EpEcalUp = new TGraph(n, x, yEup);
    TGraph *EpEcalDn = new TGraph(n, x, yEdn);
    TGraph *EpEcalL = new TGraph(n, x, yE);
    EpEcalUp->SetName("EpEcalUpperCut");
    EpEcalUp->SetLineColor(kRed);
    EpEcalDn->SetName("EpEcalLowerCut");
    EpEcalDn->SetLineColor(kRed);
    EpEcalL->SetName("EpEcalLine");
    EpEcalL->SetLineColor(kBlack);
    

    TGraph *EpHcalUp = new TGraph(n, x, yHup);
    TGraph *EpHcalDn = new TGraph(n, x, yHdn);
    TGraph *EpHcalL = new TGraph(n, x, yH);
    EpHcalUp->SetName("EpHcalUpperCut");
    EpHcalUp->SetLineColor(kRed);
    EpHcalDn->SetName("EpHcalLowerCut");
    EpHcalDn->SetLineColor(kRed);
    EpHcalL->SetName("EpHcalLine");
    EpHcalL->SetLineColor(kBlack);

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

    TCanvas *cEpTrueEcal = new TCanvas("cEpTrueEcal","cEpTrueEcal",1100,800);
    cEpTrueEcal->Divide(2,3);
    cEpTrueEcal->cd(1);
    muonEpTrueEcal->Draw("COLZ");
    EpEcalUp->Draw("SAME");
    EpEcalDn->Draw("SAME");
    EpEcalL->Draw("SAME");
    cEpTrueEcal->cd(2);
    pionEpTrueEcal->Draw("COLZ");
    EpEcalUp->Draw("SAME");
    EpEcalDn->Draw("SAME");
    EpEcalL->Draw("SAME");
    cEpTrueEcal->cd(3);
    kaonEpTrueEcal->Draw("COLZ");
    EpEcalUp->Draw("SAME");
    EpEcalDn->Draw("SAME");
    EpEcalL->Draw("SAME");
    cEpTrueEcal->cd(4);
    protonEpTrueEcal->Draw("COLZ");
    EpEcalUp->Draw("SAME");
    EpEcalDn->Draw("SAME");
    EpEcalL->Draw("SAME");
    cEpTrueEcal->cd(5);
    electronEpTrueEcal->Draw("COLZ");
    EpEcalUp->Draw("SAME");
    EpEcalDn->Draw("SAME");
    EpEcalL->Draw("SAME");
    cEpTrueEcal->Update();

    TCanvas *cEpRecEcal = new TCanvas("cEpRecEcal","cEpRecEcal",1100,800);
    cEpRecEcal->Divide(2,3);
    cEpRecEcal->cd(1);
    muonEpRecEcal->Draw("COLZ");
    EpEcalUp->Draw("SAME");
    EpEcalDn->Draw("SAME");
    EpEcalL->Draw("SAME");
    cEpRecEcal->cd(2);
    pionEpRecEcal->Draw("COLZ");
    EpEcalUp->Draw("SAME");
    EpEcalDn->Draw("SAME");
    EpEcalL->Draw("SAME");
    cEpRecEcal->cd(3);
    kaonEpRecEcal->Draw("COLZ");
    EpEcalUp->Draw("SAME");
    EpEcalDn->Draw("SAME");
    EpEcalL->Draw("SAME");
    cEpRecEcal->cd(4);
    protonEpRecEcal->Draw("COLZ");
    EpEcalUp->Draw("SAME");
    EpEcalDn->Draw("SAME");
    EpEcalL->Draw("SAME");
    cEpRecEcal->cd(5);
    electronEpRecEcal->Draw("COLZ");
    EpEcalUp->Draw("SAME");
    EpEcalDn->Draw("SAME");
    EpEcalL->Draw("SAME");
    cEpRecEcal->Update();

    TH1D *muonEpRecEcal_0 = new TH1D("muonEpRecEcal_0","Muons E/p (reco) Slices Constant;Momentum (GeV/c);Constant E/p (reco)",50,0.,25.);
    TH1D *muonEpRecEcal_1 = new TH1D("muonEpRecEcal_1","Muons E/p (reco) Slice Means;Momentum (GeV/c);E/p (reco) Mean",50,0.,25.);
    TH1D *muonEpRecEcal_2 = new TH1D("muonEpRecEcal_2","Muons E/p (reco) Slice Sigmas;Momentum (GeV/c);E/p (reco) Sigma",50,0.,25.);
    TH1D *pionEpRecEcal_0 = new TH1D("pionEpRecEcal_0","Pions E/p (reco) Slices Constant;Momentum (GeV/c);Constant E/p (reco)",50,0.,25.);
    TH1D *pionEpRecEcal_1 = new TH1D("pionEpRecEcal_1","Pions E/p (reco) Slice Means;Momentum (GeV/c);E/p (reco) Mean",50,0.,25.);
    TH1D *pionEpRecEcal_2 = new TH1D("pionEpRecEcal_2","Pions E/p (reco) Slice Sigmas;Momentum (GeV/c);E/p (reco) Sigma",50,0.,25.);

    if (muonEpRecEcal->GetEntries() > 10 &&  pionEpRecEcal->GetEntries() > 10)
    {
        piEcalSlices = true;
        muEcalSlices = true;
        muonEpRecEcal->FitSlicesY(0,0,-1,0);;
        muonEpRecEcal_0 = (TH1D*)gDirectory->Get("muonEpRecEcal_0");
        muonEpRecEcal_1 = (TH1D*)gDirectory->Get("muonEpRecEcal_1");
        muonEpRecEcal_2 = (TH1D*)gDirectory->Get("muonEpRecEcal_2");

        TCanvas *slicesMuonEpRecEcal = new TCanvas("slicesMuonEpRecEcal","slicesMuonEpRecEcal",600,600);
        slicesMuonEpRecEcal->Divide(2,2);
        slicesMuonEpRecEcal->cd(1);
        muonEpRecEcal->Draw("COLZ");
        slicesMuonEpRecEcal->cd(2);
        muonEpRecEcal_0->SetTitle("Muons E/p (reco) Slices Constant;Momentum (GeV/c);Constant E/p (reco)");
        muonEpRecEcal_0->SetLineColor(kBlue);
        muonEpRecEcal_0->Draw();
        slicesMuonEpRecEcal->cd(3);
        muonEpRecEcal_1->SetTitle("Muons E/p (reco) Slice Means;Momentum (GeV/c);E/p (reco) Mean");
        muonEpRecEcal_1->SetLineColor(kBlue);
        muonEpRecEcal_1->Draw();
        slicesMuonEpRecEcal->cd(4);
        muonEpRecEcal_2->SetTitle("Muons E/p (reco) Slice Sigmas;Momentum (GeV/c);E/p (reco) Sigma");
        muonEpRecEcal_2->SetLineColor(kBlue);
        muonEpRecEcal_2->Draw();
        slicesMuonEpRecEcal->Update();

        pionEpRecEcal->FitSlicesY(0,0,-1,0);;
        pionEpRecEcal_0 = (TH1D*)gDirectory->Get("pionEpRecEcal_0");
        pionEpRecEcal_1 = (TH1D*)gDirectory->Get("pionEpRecEcal_1");
        pionEpRecEcal_2 = (TH1D*)gDirectory->Get("pionEpRecEcal_2");

        TCanvas *slicesPionEpRecEcal = new TCanvas("slicesPionEpRecEcal","slicesPionEpRecEcal",600,600);
        slicesPionEpRecEcal->Divide(2,2);
        slicesPionEpRecEcal->cd(1);
        muonEpRecEcal->Draw("COLZ");
        slicesPionEpRecEcal->cd(2);
        pionEpRecEcal_0->SetTitle("Pions E/p (reco) Slices Constant;Momentum (GeV/c);Constant E/p (reco)");
        pionEpRecEcal_0->SetLineColor(kBlue);
        pionEpRecEcal_0->Draw();
        slicesPionEpRecEcal->cd(3);
        pionEpRecEcal_1->SetTitle("Pions E/p (reco) Slice Means;Momentum (GeV/c);E/p (reco) Mean");
        pionEpRecEcal_1->SetLineColor(kBlue);
        pionEpRecEcal_1->Draw();
        slicesPionEpRecEcal->cd(4);
        pionEpRecEcal_2->SetTitle("Pions E/p (reco) Slice Sigmas;Momentum (GeV/c);E/p (reco) Sigma");
        pionEpRecEcal_2->SetLineColor(kBlue);
        pionEpRecEcal_2->Draw();
        slicesPionEpRecEcal->Update();

    }
    else if (muonEpRecEcal->GetEntries() > pionEpRecEcal->GetEntries())
    {
        muEcalSlices = true;
        muonEpRecEcal->FitSlicesY(0,0,-1,0);;
        muonEpRecEcal_0 = (TH1D*)gDirectory->Get("muonEpRecEcal_0");
        muonEpRecEcal_1 = (TH1D*)gDirectory->Get("muonEpRecEcal_1");
        muonEpRecEcal_2 = (TH1D*)gDirectory->Get("muonEpRecEcal_2");

        TCanvas *slicesEpRecEcal = new TCanvas("slicesEpRecEcal","slicesEpRecEcal",600,600);
        slicesEpRecEcal->Divide(2,2);
        slicesEpRecEcal->cd(1);
        muonEpRecEcal->Draw("COLZ");
        slicesEpRecEcal->cd(2);
        muonEpRecEcal_0->SetTitle("Muons E/p (reco) Slices Constant;Momentum (GeV/c);Constant E/p (reco)");
        muonEpRecEcal_0->SetLineColor(kBlue);
        muonEpRecEcal_0->Draw();
        slicesEpRecEcal->cd(3);
        muonEpRecEcal_1->SetTitle("Muons E/p (reco) Slice Means;Momentum (GeV/c);E/p (reco) Mean");
        muonEpRecEcal_1->SetLineColor(kBlue);
        muonEpRecEcal_1->Draw();
        slicesEpRecEcal->cd(4);
        muonEpRecEcal_2->SetTitle("Muons E/p (reco) Slice Sigmas;Momentum (GeV/c);E/p (reco) Sigma");
        muonEpRecEcal_2->SetLineColor(kBlue);
        muonEpRecEcal_2->Draw();
        slicesEpRecEcal->Update();

    }
    else if (muonEpRecEcal->GetEntries() < pionEpRecEcal->GetEntries())
    {
        piEcalSlices = true;
        pionEpRecEcal->FitSlicesY(0,0,-1,0);;
        pionEpRecEcal_0 = (TH1D*)gDirectory->Get("pionEpRecEcal_0");
        pionEpRecEcal_1 = (TH1D*)gDirectory->Get("pionEpRecEcal_1");
        pionEpRecEcal_2 = (TH1D*)gDirectory->Get("pionEpRecEcal_2");

        TCanvas *slicesEpRecEcal = new TCanvas("slicesEpRecEcal","slicesEpRecEcal",600,600);
        slicesEpRecEcal->Divide(2,2);
        slicesEpRecEcal->cd(1);
        muonEpRecEcal->Draw("COLZ");
        slicesEpRecEcal->cd(2);
        pionEpRecEcal_0->SetTitle("Pions E/p (reco) Slices Constant;Momentum (GeV/c);Constant E/p (reco)");
        pionEpRecEcal_0->SetLineColor(kBlue);
        pionEpRecEcal_0->Draw();
        slicesEpRecEcal->cd(3);
        pionEpRecEcal_1->SetTitle("Pions E/p (reco) Slice Means;Momentum (GeV/c);E/p (reco) Mean");
        pionEpRecEcal_1->SetLineColor(kBlue);
        pionEpRecEcal_1->Draw();
        slicesEpRecEcal->cd(4);
        pionEpRecEcal_2->SetTitle("Pions E/p (reco) Slice Sigmas;Momentum (GeV/c);E/p (reco) Sigma");
        pionEpRecEcal_2->SetLineColor(kBlue);
        pionEpRecEcal_2->Draw();
        slicesEpRecEcal->Update();
    }
    


    TCanvas *cEpTrueHcal = new TCanvas("cEpTrueHcal","cEpTrueHcal",1100,800);
    cEpTrueHcal->Divide(2,3);
    cEpTrueHcal->cd(1);
    muonEpTrueHcal->Draw("COLZ");
    EpHcalUp->Draw("SAME");
    EpHcalDn->Draw("SAME");
    EpHcalL->Draw("SAME");
    cEpTrueHcal->cd(2);
    pionEpTrueHcal->Draw("COLZ");
    EpHcalUp->Draw("SAME");
    EpHcalDn->Draw("SAME");
    EpHcalL->Draw("SAME");
    cEpTrueHcal->cd(3);
    kaonEpTrueHcal->Draw("COLZ");
    EpHcalUp->Draw("SAME");
    EpHcalDn->Draw("SAME");
    EpHcalL->Draw("SAME");
    cEpTrueHcal->cd(4);
    protonEpTrueHcal->Draw("COLZ");
    EpHcalUp->Draw("SAME");
    EpHcalDn->Draw("SAME");
    EpHcalL->Draw("SAME");
    cEpTrueHcal->cd(5);
    electronEpTrueHcal->Draw("COLZ");
    EpHcalUp->Draw("SAME");
    EpHcalDn->Draw("SAME");
    EpHcalL->Draw("SAME");
    cEpTrueHcal->Update();

    TCanvas *cEpRecHcal = new TCanvas("cEpRecHcal","cEpRecHcal",1100,800);
    cEpRecHcal->Divide(2,3);
    cEpRecHcal->cd(1);
    muonEpRecHcal->Draw("COLZ");
    EpHcalUp->Draw("SAME");
    EpHcalDn->Draw("SAME");
    EpHcalL->Draw("SAME");
    cEpRecHcal->cd(2);
    pionEpRecHcal->Draw("COLZ");
    EpHcalUp->Draw("SAME");
    EpHcalDn->Draw("SAME");
    EpHcalL->Draw("SAME");
    cEpRecHcal->cd(3);
    kaonEpRecHcal->Draw("COLZ");
    EpHcalUp->Draw("SAME");
    EpHcalDn->Draw("SAME");
    EpHcalL->Draw("SAME");
    cEpRecHcal->cd(4);
    protonEpRecHcal->Draw("COLZ");
    EpHcalUp->Draw("SAME");
    EpHcalDn->Draw("SAME");
    EpHcalL->Draw("SAME");
    cEpRecHcal->cd(5);
    electronEpRecHcal->Draw("COLZ");
    EpHcalUp->Draw("SAME");
    EpHcalDn->Draw("SAME");
    EpHcalL->Draw("SAME");
    cEpRecHcal->Update();

    TH1D *muonEpRecHcal_0;
    TH1D *muonEpRecHcal_1;
    TH1D *muonEpRecHcal_2;
    TH1D *pionEpRecHcal_0;
    TH1D *pionEpRecHcal_1;
    TH1D *pionEpRecHcal_2;

    if (muonEpRecEcal->GetEntries() > 1000 &&  pionEpRecEcal->GetEntries() > 1000)
    {
        piHcalSlices = true;
        muHcalSlices = true;
        muonEpRecHcal->FitSlicesY(0,0,-1,0);;
        muonEpRecHcal_0 = (TH1D*)gDirectory->Get("muonEpRecHcal_0");
        muonEpRecHcal_1 = (TH1D*)gDirectory->Get("muonEpRecHcal_1");
        muonEpRecHcal_2 = (TH1D*)gDirectory->Get("muonEpRecHcal_2");

        TCanvas *slicesMuonEpRecHcal = new TCanvas("slicesMuonEpRecHcal","slicesMuonEpRecHcal",600,600);
        slicesMuonEpRecHcal->Divide(2,2);
        slicesMuonEpRecHcal->cd(1);
        muonEpRecHcal->Draw("COLZ");
        slicesMuonEpRecHcal->cd(2);
        muonEpRecHcal_0->SetTitle("Muons E/p (reco) Slices Constant;Momentum (GeV/c);Constant E/p (reco)");
        muonEpRecHcal_0->SetLineColor(kBlue);
        muonEpRecHcal_0->Draw();
        slicesMuonEpRecHcal->cd(3);
        muonEpRecHcal_1->SetTitle("Muons E/p (reco) Slice Means;Momentum (GeV/c);E/p (reco) Mean");
        muonEpRecHcal_1->SetLineColor(kBlue);
        muonEpRecHcal_1->Draw();
        slicesMuonEpRecHcal->cd(4);
        muonEpRecHcal_2->SetTitle("Muons E/p (reco) Slice Sigmas;Momentum (GeV/c);E/p (reco) Sigma");
        muonEpRecHcal_2->SetLineColor(kBlue);
        muonEpRecHcal_2->Draw();
        slicesMuonEpRecHcal->Update();

        pionEpRecHcal->FitSlicesY(0,0,-1,0);;
        pionEpRecHcal_0 = (TH1D*)gDirectory->Get("pionEpRecHcal_0");
        pionEpRecHcal_1 = (TH1D*)gDirectory->Get("pionEpRecHcal_1");
        pionEpRecHcal_2 = (TH1D*)gDirectory->Get("pionEpRecEcal_2");

        TCanvas *slicesPionEpRecHcal = new TCanvas("slicesPionEpRecHcal","slicesPionEpRecHcal",600,600);
        slicesPionEpRecHcal->Divide(2,2);
        slicesPionEpRecHcal->cd(1);
        muonEpRecEcal->Draw("COLZ");
        slicesPionEpRecHcal->cd(2);
        pionEpRecHcal_0->SetTitle("Pions E/p (reco) Slices Constant;Momentum (GeV/c);Constant E/p (reco)");
        pionEpRecHcal_0->SetLineColor(kBlue);
        pionEpRecHcal_0->Draw();
        slicesPionEpRecHcal->cd(3);
        pionEpRecHcal_1->SetTitle("Pions E/p (reco) Slice Means;Momentum (GeV/c);E/p (reco) Mean");
        pionEpRecHcal_1->SetLineColor(kBlue);
        pionEpRecHcal_1->Draw();
        slicesPionEpRecHcal->cd(4);
        pionEpRecHcal_2->SetTitle("Pions E/p (reco) Slice Sigmas;Momentum (GeV/c);E/p (reco) Sigma");
        pionEpRecHcal_2->SetLineColor(kBlue);
        pionEpRecHcal_2->Draw();
        slicesPionEpRecHcal->Update();

    }
    else if (muonEpRecEcal->GetEntries() > pionEpRecEcal->GetEntries())
    {
        muHcalSlices = true;
        muonEpRecHcal->FitSlicesY(0,0,-1,0);;
        muonEpRecHcal_0 = (TH1D*)gDirectory->Get("muonEpRecHcal_0");
        muonEpRecHcal_1 = (TH1D*)gDirectory->Get("muonEpRecHcal_1");
        muonEpRecHcal_2 = (TH1D*)gDirectory->Get("muonEpRecHcal_2");

        TCanvas *slicesMuonEpRecHcal = new TCanvas("slicesMuonEpRecHcal","slicesMuonEpRecHcal",600,600);
        slicesMuonEpRecHcal->Divide(2,2);
        slicesMuonEpRecHcal->cd(1);
        muonEpRecHcal->Draw("COLZ");
        slicesMuonEpRecHcal->cd(2);
        muonEpRecHcal_0->SetTitle("Muons E/p (reco) Slices Constant;Momentum (GeV/c);Constant E/p (reco)");
        muonEpRecHcal_0->SetLineColor(kBlue);
        muonEpRecHcal_0->Draw();
        slicesMuonEpRecHcal->cd(3);
        muonEpRecHcal_1->SetTitle("Muons E/p (reco) Slice Means;Momentum (GeV/c);E/p (reco) Mean");
        muonEpRecHcal_1->SetLineColor(kBlue);
        muonEpRecHcal_1->Draw();
        slicesMuonEpRecHcal->cd(4);
        muonEpRecHcal_2->SetTitle("Muons E/p (reco) Slice Sigmas;Momentum (GeV/c);E/p (reco) Sigma");
        muonEpRecHcal_2->SetLineColor(kBlue);
        muonEpRecHcal_2->Draw();
        slicesMuonEpRecHcal->Update();

    }
    else if (muonEpRecEcal->GetEntries() < pionEpRecEcal->GetEntries())
    {
        piHcalSlices = true;
        pionEpRecHcal->FitSlicesY(0,0,-1,0);;
        pionEpRecHcal_0 = (TH1D*)gDirectory->Get("pionEpRecHcal_0");
        pionEpRecHcal_1 = (TH1D*)gDirectory->Get("pionEpRecHcal_1");
        pionEpRecHcal_2 = (TH1D*)gDirectory->Get("pionEpRecEcal_2");

        TCanvas *slicesPionEpRecHcal = new TCanvas("slicesPionEpRecHcal","slicesPionEpRecHcal",600,600);
        slicesPionEpRecHcal->Divide(2,2);
        slicesPionEpRecHcal->cd(1);
        muonEpRecEcal->Draw("COLZ");
        slicesPionEpRecHcal->cd(2);
        pionEpRecHcal_0->SetTitle("Pions E/p (reco) Slices Constant;Momentum (GeV/c);Constant E/p (reco)");
        pionEpRecHcal_0->SetLineColor(kBlue);
        pionEpRecHcal_0->Draw();
        slicesPionEpRecHcal->cd(3);
        pionEpRecHcal_1->SetTitle("Pions E/p (reco) Slice Means;Momentum (GeV/c);E/p (reco) Mean");
        pionEpRecHcal_1->SetLineColor(kBlue);
        pionEpRecHcal_1->Draw();
        slicesPionEpRecHcal->cd(4);
        pionEpRecHcal_2->SetTitle("Pions E/p (reco) Slice Sigmas;Momentum (GeV/c);E/p (reco) Sigma");
        pionEpRecHcal_2->SetLineColor(kBlue);
        pionEpRecHcal_2->Draw();
        slicesPionEpRecHcal->Update();
    }

    TCanvas *cEpVsEta = new TCanvas("cEpVsEta","cEpVsEta",1100,800);
    cEpVsEta->Divide(2,2);
    cEpVsEta->cd(1);;
    mupiEpTrueEcalvsEta->Draw("COLZ");
    cEpVsEta->cd(2);
    mupiEpRecEcalvsEta->Draw("COLZ");
    cEpVsEta->cd(3);
    mupiEpTrueHcalvsEta->Draw("COLZ");
    cEpVsEta->cd(4);
    mupiEpRecHcalvsEta->Draw("COLZ");
    cEpVsEta->Update();

    TCanvas *cEpTrueEHcal = new TCanvas("cEpTrueEHcal","cEpTrueEHcal",1100,800);
    cEpTrueEHcal->Divide(2,3);
    cEpTrueEHcal->cd(1);
    muonEpTrueEHcal->Draw("COLZ");
    cEpTrueEHcal->cd(2);
    pionEpTrueEHcal->Draw("COLZ");
    cEpTrueEHcal->cd(3);
    kaonEpTrueEHcal->Draw("COLZ");
    cEpTrueEHcal->cd(4);
    protonEpTrueEHcal->Draw("COLZ");
    cEpTrueEHcal->cd(5);
    electronEpTrueEHcal->Draw("COLZ");
    cEpTrueEHcal->Update();

    TCanvas *cEpRecEHcal = new TCanvas("cEpRecEHcal","cEpRecEHcal",1100,800);
    cEpRecEHcal->Divide(2,3);
    cEpRecEHcal->cd(1);
    muonEpRecEHcal->Draw("COLZ");
    cEpRecEHcal->cd(2);
    pionEpRecEHcal->Draw("COLZ");
    cEpRecEHcal->cd(3);
    kaonEpRecEHcal->Draw("COLZ");
    cEpRecEHcal->cd(4);
    protonEpRecEHcal->Draw("COLZ");
    cEpRecEHcal->cd(5);
    electronEpRecEHcal->Draw("COLZ");
    cEpRecEHcal->Update();

    TCanvas *cEcalClusterSize = new TCanvas("cEcalClusterSize","cEcalClusterSize",1100,800);
    cEcalClusterSize->Divide(2,3);
    cEcalClusterSize->cd(1);
    muonEcalClusterSize->Draw("colz");
    cEcalClusterSize->cd(2);
    pionEcalClusterSize->Draw("colz");
    cEcalClusterSize->cd(3);
    kaonEcalClusterSize->Draw("colz");
    cEcalClusterSize->cd(4);
    protonEcalClusterSize->Draw("colz");
    cEcalClusterSize->cd(5);
    electronEcalClusterSize->Draw("colz");
    cEcalClusterSize->Update();

    TCanvas *cHcalClusterSize = new TCanvas("cHcalClusterSize","cHcalClusterSize",1100,800);
    cHcalClusterSize->Divide(2,3);
    cHcalClusterSize->cd(1);
    muonHcalClusterSize->Draw("colz");
    cHcalClusterSize->cd(2);
    pionHcalClusterSize->Draw("colz");
    cHcalClusterSize->cd(3);
    kaonHcalClusterSize->Draw("colz");
    cHcalClusterSize->cd(4);
    protonHcalClusterSize->Draw("colz");
    cHcalClusterSize->cd(5);
    electronHcalClusterSize->Draw("colz");
    cHcalClusterSize->Update();

    TCanvas *cEHcalClusterSize = new TCanvas("cEHcalClusterSize","cEHcalClusterSize",1100,800);
    cEHcalClusterSize->Divide(2,3);
    cEHcalClusterSize->cd(1);
    muonEHcalClusterSize->Draw("colz");
    cEHcalClusterSize->cd(2);
    pionEHcalClusterSize->Draw("colz");
    cEHcalClusterSize->cd(3);
    kaonEHcalClusterSize->Draw("colz");
    cEHcalClusterSize->cd(4);
    protonEHcalClusterSize->Draw("colz");
    cEHcalClusterSize->cd(5);
    electronEHcalClusterSize->Draw("colz");
    cEHcalClusterSize->Update();


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


    muonMomEff->Divide(muonTrueP);
    muonMomPur->Divide(allParticlesTrueP_PC);

    muonMomEtaEff->Divide(muonTruePEta);
    muonMomEtaPur->Divide(allParticlesTrueEtaP_PC);

    TCanvas *cMuonEffPur = new TCanvas("cMuonEffPur","cMuonEffPur",1100,800);
    cMuonEffPur->Divide(2,2);
    cMuonEffPur->cd(1);
    muonMomEff->Draw("COLZ");
    cMuonEffPur->cd(2);
    muonMomPur->Draw("COLZ");
    cMuonEffPur->cd(3);
    muonMomEtaEff->Draw("COLZ");
    cMuonEffPur->cd(4);
    muonMomEtaPur->Draw("COLZ");
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
    cPostCuts->cd(5);
    muonRecMinusTrueP->Draw();
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
    muonTrueP->Write();
    ofile->cd("..");
    ofile->mkdir("trueEcalEp");
    ofile->cd("trueEcalEp");
    EpEcalUp->Write();
    EpEcalDn->Write();
    EpEcalL->Write();
    protonEpTrueEcal->Write();
    electronEpTrueEcal->Write();
    pionEpTrueEcal->Write();
    kaonEpTrueEcal->Write();
    muonEpTrueEcal->Write();
    ofile->cd("..");
    ofile->mkdir("recoEcalEp");
    ofile->cd("recoEcalEp");
    EpEcalUp->Write();
    EpEcalDn->Write();
    EpEcalL->Write();
    protonEpRecEcal->Write();
    electronEpRecEcal->Write();
    pionEpRecEcal->Write();
    if (pionEpRecEcal_0 != nullptr)
    {
        pionEpRecEcal_0->Write();
        pionEpRecEcal_1->Write();
        pionEpRecEcal_2->Write();
    }
    kaonEpRecEcal->Write();
    muonEpRecEcal->Write();
    if (muonEpRecEcal_0 != nullptr)
    {
        muonEpRecEcal_0->Write();
        muonEpRecEcal_1->Write();
        muonEpRecEcal_2->Write();
    }
    ofile->cd("..");
    ofile->mkdir("trueHcalEp");
    ofile->cd("trueHcalEp");
    EpHcalUp->Write();
    EpHcalDn->Write();
    EpHcalL->Write();
    protonEpTrueHcal->Write();
    electronEpTrueHcal->Write();
    pionEpTrueHcal->Write();
    kaonEpTrueHcal->Write();
    muonEpTrueHcal->Write();
    ofile->cd("..");
    ofile->mkdir("recoHcalEp");
    ofile->cd("recoHcalEp");
    EpHcalUp->Write();
    EpHcalDn->Write();
    EpHcalL->Write();
    protonEpRecHcal->Write();
    electronEpRecHcal->Write();
    pionEpRecHcal->Write();
    if (pionEpRecHcal_0 != nullptr)
    {
        pionEpRecHcal_0->Write();
        pionEpRecHcal_1->Write();
        pionEpRecHcal_2->Write();
    }
    kaonEpRecHcal->Write();
    muonEpRecHcal->Write();
    if (muonEpRecHcal_0 != nullptr)
    {
        muonEpRecHcal_0->Write();
        muonEpRecHcal_1->Write();
        muonEpRecHcal_2->Write();
    }
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
    muonRecMinusTrueP->Write();
    pionPostCutEta->Write();
    pionPostCutP->Write();
    muonMomEff->Write();
    muonMomEtaEff->Write();
    muonMomPur->Write();
    muonMomEtaPur->Write();
    ofile->cd("..");


}