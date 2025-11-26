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
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLatex.h>
#include <Math/Boost.h>

#include "ePICStyle.C"
#include "MuonFinder.C"
#include "ElectronFinder.C"
#include "ParentFinder.C"
#include "DVMP_JPsi_Analysis.h"

void DVMP_JPsi_Analysis()
{
    gROOT->SetBatch(kTRUE);
    gROOT->ProcessLine("SetePICStyle()");
    //gStyle->SetOptStat(0);

    TString infile="eicReconOutput/EICreconOut_JPsiMuMu_10ifb_10x250ep_Pruned.root";
    //TString infile="dis_background/DIS_Q2_1_10_10x130ep_Pruned.root";
    
    std::string filename = infile.Data();
    std::string beam_config;
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

    double protonEnergy = std::stod(proton_energy); 
    double lumi = 10.0; // Luminosity in fb^-1, adjust as needed

    std::cout << "Extracted beam config: " << beam_config << std::endl;
    std::cout << "Proton energy: " << protonEnergy << std::endl;

    int lumi_int = static_cast<int>(lumi);
    std::string outfilename = "outputs/DVMP_JPsi_AnalysisOutput_" + std::to_string(lumi_int) + "ifb_" + beam_config + ".root";
    //std::string outfilename = "outputs/DIS_Q2_1_10_AnalysisOutput_" + beam_config + ".root";

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

    // Get Forward Detectors Information
    TTreeReaderArray<float> RPEng(tree_reader, "ForwardRomanPotRecParticles.energy");
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

    // Define Eta Histograms

    TH1D *electronEta = new TH1D("electronEta", " #eta of Thrown Electrons",150,-15.,15.);
    TH1D *matchedElectronEta = new TH1D("matchedElectronEta", " #eta of Thrown Electrons That Have Matching Track",150,-15.,15.);
    TH1D *electronEff = new TH1D("electronEff","Efficency;  #eta",150,-15.,15.);

    TH1D *protonEta = new TH1D("protonEta", " #eta of Thrown Protons",150,-15.,15.);
    TH1D *matchedProtonEta = new TH1D("matchedProtonEta", " #eta of Thrown Protons That Have Matching Track",150,-15.,15.);
    TH1D *protonEff = new TH1D("protonEff","Efficency;  #eta",150,-15.,15.);

    TH1D *muonEta = new TH1D("muonEta", " #eta of Thrown Muons",100,-5.,5.);
    TH1D *matchedMuonEta = new TH1D("matchedMuonEta", " #eta of Thrown Muons That Have Matching Track",100,-5.,5.);
    TH1D *muonEff = new TH1D("muonEff","Efficency;  #eta",100,-5.,5.);

    TH1D *JPsiEta = new TH1D("JPsiEta", " #eta of Thrown J/PSi",100,-5.,5.);
    TH1D *matchedJPsiEta = new TH1D("matchedJPsiEta", " #eta of Thrown J/PSi That Have Matching Track",100,-5.,5.);
    TH1D *JPsiEff = new TH1D("JPsiEff","Efficency;  #eta",100,-5.,5.);

    // Define Momentum Histograms

    TH1D *electronMomHist = new TH1D("electronMomHist","Momentum of Thrown Electrons;Momentum (GeV/c)",100,-100.,100.);
    TH1D *matchedElectronMomHist = new TH1D("matchedElectronMomHist","Momentum of Thrown Electrons That Have Matching Track (GeV/c)",100,-100.,100.);
    TH1D *electronMomEff = new TH1D("electronMomEff","Efficency;Momentum (GeV/c)",100,-100.,100.);

    TH1D *protonMomHist = new TH1D("protonMomHist","Momentum of Thrown Protons;Momentum (GeV/c)",160,-20.,300);
    TH1D *matchedProtonMomHist = new TH1D("matchedProtonMomHist","Momentum of Thrown Protons That Have Matching Track (GeV/c)" ,160,-20.,300);
    TH1D *protonMomEff = new TH1D("protonMomEff","Efficency;Momentum (GeV/c)",160,-20.,300);

    TH1D *protonEMinusPzHist = new TH1D("protonEMinusPzHist","E - Pz of Thrown Protons;E - Pz (GeV)",160,-20.,300);
    TH1D *matchedProtonEMinusPzHist = new TH1D("matchedProtonEMinusPzHist","E - Pz of Thrown Protons That Have Matching Track;E - Pz (GeV)",160,-20.,300);
    TH1D *protonPtHist = new TH1D("protonPtHist","Pt of Thrown Protons;Pt (GeV/c)",100,0.,10.);
    TH1D *matchedProtonPtHist = new TH1D("matchedProtonPtHist","Pt of Thrown Protons That Have Matching Track;Pt (GeV/c)",100,0.,10.);

    TH1D *electronProtonEMinusPzHist = new TH1D("electronProtonEMinusPzHist","E - Pz of Thrown e+p System;E - Pz (GeV)",200,0.,200.);
    TH1D *matchedElectronProtonEMinusPzHist = new TH1D("matchedElectronProtonEMinusPzHist","E - Pz of Thrown e+p System That Have Matching Tracks;E - Pz (GeV)",200,0.,200.);
    TH1D *electronProtonPtHist = new TH1D("electronProtonPtHist","Pt of Thrown e+p System;Pt (GeV/c)",100,0.,10.);
    TH1D *matchedElectronProtonPtHist = new TH1D("matchedElectronProtonPtHist","Pt of Thrown e+p System That Have Matching Tracks;Pt (GeV/c)",100,0.,10.);

    TH1D *muonMomHist = new TH1D("muonMomHist","Momentum of Thrown Muons;Momentum (GeV/c)",160,-20.,300);
    TH1D *matchedMuonMomHist = new TH1D("matchedMuonMomHist","Momentum of Thrown Muons That Have Matching Track (GeV/c)",160,-20.,300);
    TH1D *muonMomEff = new TH1D("muonMomEff","Efficency;Momentum (GeV/c)",160,-20.,300);

    TH1D *JPsiMomHist = new TH1D("JPsiMomHist","Momentum of Thrown J/PSi;Momentum (GeV/c)",160,-20.,300);
    TH1D *matchedJPsiMomHist = new TH1D("matchedJPsiMomHist","Momentum of Thrown J/PSi That Have Matching Track (GeV/c)",160,-20.,300);
    TH1D *JPsiMomEff = new TH1D("JPsiMomEff","Efficency;Momentum (GeV/c)",160,-20.,300);

    // Define 2D Histograms p vs eta
    TH2D *electronEtaMom = new TH2D("electronEtaMom","Thrown Electrons;  #eta; Momentum (GeV/c)",150,-15.,15.,100,-100.,100.);
    TH2D *matchedElectronEtaMom = new TH2D("matchedElectronEtaMom","Thrown Electrons That Have Matching Track;  #eta; Momentum (GeV/c)",150,-15.,15.,100,-100.,100.);

    TH2D *protonEtaMom = new TH2D("protonEtaMom","Thrown Protons;  #eta; Momentum (GeV/c)",150,-15.,15.,160,-20.,300);
    TH2D *matchedProtonEtaMom = new TH2D("matchedProtonEtaMom","Thrown Protons That Have Matching Track;  #eta; Momentum (GeV/c)",150,-15.,15.,160,-20.,300);

    TH2D *muonEtaMom = new TH2D("muonEtaMom","Thrown Muons;  #eta; Momentum (GeV/c)",200,-10.,10.,160,-20.,300);
    TH2D *matchedMuonEtaMom = new TH2D("matchedMuonEtaMom","Thrown Muons That Have Matching Track;  #eta; Momentum (GeV/c)",200,-10.,10.,160,-20.,300);

    TH2D *JPsiEtaMom = new TH2D("JPsiEtaMom","Thrown J/PSi;  #eta; Momentum (GeV/c)",100,-5.,5.,160,-20.,300);
    TH2D *matchedJPsiEtaMom = new TH2D("matchedJPsiEtaMom","Thrown J/PSi That Have Matching Track;  #eta; Momentum (GeV/c)",100,-5.,5.,160,-20.,300);

    // Define Delta R Histograms

    TH1D *matchedElectronTrackDeltaR = new TH1D("matchedElectronTrackDeltaR","Delta R Between Matching Thrown and Reconstructed Electron",5000,0.,5.);
    TH1D *matchedProtonTrackDeltaR = new TH1D("matchedProtonTrackDeltaR","Delta R Between Matching Thrown and Reconstructed Proton",5000,0.,5.);
    TH1D *matchedMuonTrackDeltaR = new TH1D("matchedMuonTrackDeltaR","Delta R Between Matching Thrown and Reconstructed Muon",5000,0.,5.);
    TH1D *matchedJPsiTrackDeltaR = new TH1D("matchedJPsiTrackDeltaR","Delta R Between Matching Thrown and Reconstructed J/PSi",5000,0.,5.);

    // Define Invariant Mass, Missing Mass, and Angle Histograms

    TH1D *muonScatterAngle = new TH1D("muonScatterAngle","Scattering Angle of Thrown Muons; Scattering Angle (rad)",100,0.,7.);

    TH1D *invMassMuMu = new TH1D("invMassMuMu","Invariant Mass of MC Muon Pairs; Invariant Mass (GeV)",100,0.,5.); // Invariant mass histogram
    TH1D *matchedMuMuInvMass = new TH1D("matchedMuMuInvMass","Invariant Mass of Reconstructed Muon Pairs; Invariant Mass (GeV)",100,0.,5.); // Invariant mass histogram

    TH1D *matchedMissingMassEpToJPsiX = new TH1D("matchedMissingMassEpToJPsiX","Missing Mass in e+p -> e' + J/Psi + p' from matched tracks; Missing Mass (GeV/c^{2})",400,-20.,20.);
    TH1D *matchedMissingMass2EpToJPsiX = new TH1D("matchedMissingMass2EpToJPsiX","Missing Mass Squared in e+p -> e' + J/Psi + p' from matched tracks; Missing Mass Squared (GeV/c^{2})^{2}",800,-40.,40.);

    // Define Kinematics Histograms
    TH1D* trueQ2 = new TH1D("trueQ^{2}","True Q^{2} Distribution; Q^{2} (GeV/c^{2})",100,0.,50.);
    TH1D* reconQ2_DA = new TH1D("reconQ^{2}_DA","Reconstructed Q^{2} Distribution using the DA Method; Q^{2} (GeV/c^{2})",100,0.,50.);
    TH1D* reconQ2_JB = new TH1D("reconQ^{2}_JB","Reconstructed Q^{2} Distribution using the JB Method; Q^{2} (GeV/c^{2})",100,0.,50.);
    TH1D* reconQ2_e = new TH1D("reconQ^{2}_e","Reconstructed Q^{2} Distribution using the electron Method; Q^{2} (GeV/c^{2})",100,0.,50.);
    TH1D* reconQ2_sigma = new TH1D("reconQ^{2}_sigma","Reconstructed Q^{2} Distribution using the sigma Method; Q^{2} (GeV/c^{2})",100,0.,50.);

    TH1D* deltaQ2_DA = new TH1D("deltaQ2_DA","Delta Q^{2} (Reconstructed - True) using the DA Method; (Q^{2} - Q_{MC}^{2})/Q_{MC}^{2} (%)",200,-100.,100.);
    TH1D* deltaQ2_JB = new TH1D("deltaQ2_JB","Delta Q^{2} (Reconstructed - True) using the JB Method; (Q^{2} - Q_{MC}^{2})/Q_{MC}^{2} (%)",200,-100.,100.);
    TH1D* deltaQ2_e = new TH1D("deltaQ2_e","Delta Q^{2} (Reconstructed - True) using the electron Method; (Q^{2} - Q_{MC}^{2})/Q_{MC}^{2} (%)",200,-100.,100.);
    TH1D* deltaQ2_sigma = new TH1D("deltaQ2_sigma","Delta Q^{2} (Reconstructed - True) using the sigma Method; (Q^{2} - Q_{MC}^{2})/Q_{MC}^{2} (%)",200,-100.,100.);

    TH2D* reconQ2_DA_vs_trueQ2 = new TH2D("reconQ2_DA_vs_trueQ2","Reconstructed Q^{2} using the DA Method vs True Q^{2}; Q_{MC}^{2} (GeV/c^{2}); Q_{DA}^{2} (GeV/c^{2})",100,0.,50.,100,0.,50.);

    TH1D* truet = new TH1D("truet","True t Distribution; t (GeV/c)",100,0.,2.);
    TH1D* recont_eXBABE = new TH1D("recont_eXBABE","Reconstructed t Distribution using the eXBABE Method; t (GeV/c)",100,0.,2.);
    TH1D* recont_eXPT = new TH1D("recont_eXPT","Reconstructed t Distribution using the eXPT Method; t (GeV/c)",100,0.,2.);
    TH1D* recont_eX = new TH1D("recont_eX","Reconstructed t Distribution using the eX Method; t (GeV/c)",100,0.,2.);
    TH1D* recont_BABE = new TH1D("recont_BABE","Reconstructed t Distribution using the BABE Method; t (GeV/c)",100,0.,2.);

    TH1D* deltat_eXBABE = new TH1D("deltat_eXBABE","Delta t (Reconstructed - True) using the eXBABE Method; (t - t_{MC})/t_{MC} (%)",200,-100.,100.);
    TH1D* deltat_eXPT = new TH1D("deltat_eXPT","Delta t (Reconstructed - True) using the eXPT Method; (t - t_{MC})/t_{MC} (%)",200,-100.,100.);
    TH1D* deltat_eX = new TH1D("deltat_eX","Delta t (Reconstructed - True) using the eX Method; (t - t_{MC})/t_{MC} (%)",200,-100.,100.);
    TH1D* deltat_BABE = new TH1D("deltat_BABE","Delta t (Reconstructed - True) using the BABE Method; (t - t_{MC})/t_{MC} (%)",200,-100.,100.);

    TH2D* recont_eXBABE_vs_truet = new TH2D("recont_eXBABE_vs_truet","Reconstructed t using the eXBABE Method vs True t; t_{MC} (GeV/c); t_{eXBABE} (GeV/c)",100,0.,2.,100,0.,2.);
    TH2D* deltat_eXBABE_vs_truet = new TH2D("deltat_eXBABE_vs_truet","Delta t (Reconstructed - True) using the eXBABE Method vs True t; t_{MC} (GeV/c); (t_{eXBABE} - t_{MC})/t_{MC}",100,0.,2.,200,-100.,100.);

    TH1D* truey = new TH1D("truey","True y Distribution",100,0.,1.0);
    TH1D* recony_DA = new TH1D("recony_DA","Reconstructed y Distribution using the DA method; y",100,0.,1.0);
    TH1D* recony_JB = new TH1D("recony_JB","Reconstructed x Distribution using the JB method; y",100,0.,1.0);
    TH1D* recony_e = new TH1D("recony_e","Reconstructed x Distribution using the electron method; y",100,0.,1.0);
    TH1D* recony_sigma = new TH1D("recony_sigma","Reconstructed x Distribution using the sigma method; y",100,0.,1.0);

    TH2D* recony_DA_vs_truey = new TH2D("recony_DA_vs_truey","Reconstructed y using the DA Method vs True y; y_{MC}; y_{DA}^{2}",100,0.,1.,100,0.,1.);

    TH1D* deltay_DA = new TH1D("deltay_DA","Delta y (Reconstructed - True) using the DA Method; (y - y_{MC})/y_{MC} (%)",200,-100.,100.);
    TH1D* deltay_JB = new TH1D("deltay_JB","Delta y (Reconstructed - True) using the JB Method; (y - y_{MC})/y_{MC} (%)",200,-100.,100.);
    TH1D* deltay_e = new TH1D("deltay_e","Delta y (Reconstructed - True) using the electron Method; (y - y_{MC})/y_{MC} (%)",200,-100.,100.);
    TH1D* deltay_sigma = new TH1D("deltay_sigma","Delta y (Reconstructed - True) using the sigma Method; (y - y_{MC})/y_{MC} (%)",200,-100.,100.);

    TH1D* truex = new TH1D("truex","True x Distribution; x_bjk",100,0.,0.25);
    TH1D* reconx_DA = new TH1D("reconx_DA","Reconstructed x Distribution using the DA method; x_bjk",100,0.,0.25);
    TH1D* reconx_JB = new TH1D("reconx_JB","Reconstructed x Distribution using the JB method; x_bjk",100,0.,0.25);
    TH1D* reconx_e = new TH1D("reconx_e","Reconstructed x Distribution using the electron method; x_bjk",100,0.,0.25);
    TH1D* reconx_sigma = new TH1D("reconx_sigma","Reconstructed x Distribution using the sigma method; x_bjk",100,0.,0.25);

    TH1D* deltax_DA = new TH1D("deltax_DA","Delta x (Reconstructed - True) using the DA Method; (x - x_{MC})/x_{MC} (%)",200,-100.,100.);
    TH1D* deltax_JB = new TH1D("deltax_JB","Delta x (Reconstructed - True) using the JB Method; (x - x_{MC})/x_{MC} (%)",200,-100.,100.);
    TH1D* deltax_e = new TH1D("deltax_e","Delta x (Reconstructed - True) using the electron Method; (x - x_{MC})/x_{MC} (%)",200,-100.,100.);
    TH1D* deltax_sigma = new TH1D("deltax_sigma","Delta x (Reconstructed - True) using the sigma Method; (x - x_{MC})/x_{MC} (%)",200,-100.,100.);

    TH1D* truet_XbjkA = new TH1D("truet_XbjkA","True t distribution with bjorken x binning; t_{MC} (GeV/c)",100,0.,2.);
    TH1D* recont_XbjkA = new TH1D("recont_XbjkA","Reconstructed t distribution with bjorken x binning; t_{MC} (GeV/c)",100,0.,2.);
    TH1D* truet_XbjkB = new TH1D("truet_XbjkB","True t distribution with bjorken x binning; t_{MC} (GeV/c)",100,0.,2.);
    TH1D* recont_XbjkB = new TH1D("recont_XbjkB","Reconstructed t distribution with bjorken x binning; t_{eXBABE} (GeV/c)",100,0.,2.);
    TH1D* truet_XbjkC = new TH1D("truet_XbjkC","True t distribution with bjorken x binning; t_{eXBABE} (GeV/c)",100,0.,2.);
    TH1D* recont_XbjkC = new TH1D("recont_XbjkC","Reconstructed t distribution with bjorken x binning; t_{eXBABE} (GeV/c)",100,0.,2.);

    
    bool eventFail = false; 
    bool electronFound = false; // Flag to check if scattered electron is found
    bool electronFinderReturn = false; // Return value from electron finder
    bool muonFound = false; // Flag to check if muon is found
    bool muonFinderReturn = false; // Return value from muon finder
    bool invMassError = false;
    int cutEvents[6] = {0,0,0,0,0,0}; // Cut flow counter
    int eventsPassed = 0; // Events passed counter
    int eventID = 0; // Event ID counter

    int electronFinderInconclusiveRate = 0;
    int electronFinderCorrectRate = 0;
    int electronFinderIncorrectRate = 0;
    int muonFinderInconclusiveRate = 0;
    int muonFinderCorrectRate = 0;
    int muonFinderIncorrectRate = 0; 
    int incorrectElectronIDrate = 0; 
    int incorrectMuonIDrate = 0;
    int incorrectParentIDrate = 0;

    while(tree_reader.Next()) // Loop over events
    {
        eventID++;

        if (eventID % 10 == 0)
        {
            fprintf (stderr, "%4.2f Percent\r ", eventID*100.0/mychain->GetEntries());
            fflush (stderr);
        }

        if (trackEng.GetSize() != 3) // Only look at events with 3 tracks
        { 
            cutEvents[0]++;
            continue; 
        }
        if ((trackEng[0] == 0. || trackEng[1] == 0. || trackEng[2] == 0.) && (trackCharge[0] + trackCharge[1] + trackCharge[2] != 1.)) // Remove events with zero energy tracks and wrong charge sum
        { 
            cutEvents[0]++;
            continue; 
        }

        // Reset vectors
        for (int i=0; i<3; i++)
        {
            recoTrackMom[i] = TVector3(0.,0.,0.);
            recoTrack4Mom[i].SetPxPyPzE(0.,0.,0.,0.);
            recoTrackIndex[i] = -1;
            recoTrackTruePID[i] = -1;
        }

        // Reset all momentum vectors
        std::vector<TVector3*> vec3D = {&scatEMomT, &scatEMomR, &scatpMomT, &scatpMomR, 
            &muPlusMomT, &muPlusMomR, &muMinusMomT, &muMinusMomR, &JPsiMomT, &JPsiMomR};
        std::vector<ROOT::Math::PxPyPzEVector*> vec4D = {&scatE4MomT, &scatE4MomR, &scatp4MomT, &scatp4MomR,
            &muPlus4MomT, &muPlus4MomR, &muMinus4MomT, &muMinus4MomR, &JPsi4MomT, &JPsi4MomR};
            
        for(auto& v : vec3D) v->SetXYZ(0,0,0);
        for(auto& v : vec4D) v->SetXYZT(0,0,0,0);

        electronFound = false;
        electronFinderReturn = false;
        muonFound = false;
        muonFinderReturn = false;
        invMassError = false;
        eventFail = false;

        // Reset kinematic variables
        {
            double* ks[] = { &Q2_truth, &Q2_DA, &Q2_JB, &Q2_e, &Q2_sigma,
                     &t_truth, &t_eXBABE, &t_eXPT, &t_eX, &t_BABE,
                     &x_truth, &x_DA, &x_JB, &x_e, &x_sigma,
                     &y_truth, &y_DA, &y_JB, &y_e, &y_sigma };
            for (double* p : ks) *p = 0.;
        }

        // Find the true particles

        for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop over thrown particles
        {
            if (partGenStat[i] == 4) // Look for beam particles
            {
                int pdg = TMath::Abs(partPdg[i]);
                if(pdg == 11) // Electron beam
                {
                    beamEMom.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);
                    double beamEEng = sqrt(pow(partMomX[i],2) + pow(partMomY[i],2) + pow(partMomZ[i],2) + pow(partMass[i],2));
                    beamE4Mom.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], beamEEng);
                }
                if(pdg == 2212) // Proton beam
                {
                    beampMom.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);
                    beampMom.RotateY(crossingAngle);
                    double beampEng = sqrt(pow(partMomX[i],2) + pow(partMomY[i],2) + pow(partMomZ[i],2) + pow(partMass[i],2));
                    beamp4Mom.SetPxPyPzE(beampMom.X(), beampMom.Y(), beampMom.Z(), beampEng);
                }
            }

            if (partGenStat[i] == 1) // Select stable thrown particles
            {
                int pdg = TMath::Abs(partPdg[i]);
                double partEng = sqrt(pow(partMomX[i],2) + pow(partMomY[i],2) + pow(partMomZ[i],2) + pow(partMass[i],2)); // Energy of Monte Carlo particles

                if(pdg == 11) // Look at electrons
                {
                    scatEMomT = TVector3(partMomX[i],partMomY[i],partMomZ[i]);
                    scatE4MomT.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);
                    electronEta->Fill(scatEMomT.PseudoRapidity());
                    electronMomHist->Fill(scatEMomT.Mag());
                    electronEtaMom->Fill(scatEMomT.PseudoRapidity(),scatEMomT.Mag());
                }
                if(pdg == 2212) // Look at protons
                {
                    scatpMomT = TVector3(partMomX[i],partMomY[i],partMomZ[i]);
                    scatpMomT.RotateY(crossingAngle); 
                    scatp4MomT.SetPxPyPzE(scatpMomT.X(),scatpMomT.Y(),scatpMomT.Z(), partEng);
                    protonEta->Fill(scatpMomT.PseudoRapidity());
                    protonMomHist->Fill(scatpMomT.Mag());
                    protonEtaMom->Fill(scatpMomT.PseudoRapidity(),scatpMomT.Mag());
                    protonEMinusPzHist->Fill(partEng - scatpMomT.Z());
                    protonPtHist->Fill(scatpMomT.Perp());
                }
                if(pdg == 13 && partCharge[i] == 1) // Look at mu+
                {
                    muPlusMomT = TVector3(partMomX[i],partMomY[i],partMomZ[i]);
                    muPlus4MomT.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);
                    muonEta->Fill(muPlusMomT.PseudoRapidity());
                    muonMomHist->Fill(muPlusMomT.Mag());
                    muonEtaMom->Fill(muPlusMomT.PseudoRapidity(),muPlusMomT.Mag());
                }
                if(pdg == 13 && partCharge[i] == -1) // Look at mu-
                {
                    muMinusMomT = TVector3(partMomX[i],partMomY[i],partMomZ[i]);
                    muMinus4MomT.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);
                    muonEta->Fill(muMinusMomT.PseudoRapidity());
                    muonMomHist->Fill(muMinusMomT.Mag());
                    muonEtaMom->Fill(muMinusMomT.PseudoRapidity(),muMinusMomT.Mag());
                }
            }
        }
        
        if (scatEMomT.Mag() != 0. && scatpMomT.Mag() != 0.)
        { 
            electronProtonEMinusPzHist->Fill(scatE4MomT.E() + scatp4MomT.E() - (scatE4MomT.Pz() + scatp4MomT.Pz()));
            electronProtonPtHist->Fill((scatEMomT + scatpMomT).Perp());
        }

        if (muPlusMomT.Mag() != 0. && muMinusMomT.Mag() != 0.) 
        { 
            JPsiMomT = muPlusMomT + muMinusMomT; 
            JPsi4MomT = muPlus4MomT + muMinus4MomT; 

            invMassMuMu->Fill((muPlus4MomT + muMinus4MomT).M());

            JPsiEta->Fill(JPsiMomT.PseudoRapidity());
            JPsiMomHist->Fill(JPsiMomT.Mag());
            JPsiEtaMom->Fill(JPsiMomT.PseudoRapidity(),JPsiMomT.Mag());
        }
        

        // Search for the proton in the forward detectors
        for (unsigned int i = 0; i < RPEng.GetSize(); i++)
        {
            if (RPEng[i] < protonEnergy*0.5) continue; // Only look at protons with at least 50% of the beam energy
            if (RPMomZ[i] < 0.) RPMomZ[i] = -RPMomZ[i]; // Correct for negative z momentum in RP
            scatpMomR = TVector3(RPMomX[i],RPMomY[i],RPMomZ[i]);
            scatp4MomR.SetPxPyPzE(RPMomX[i],RPMomY[i],RPMomZ[i], RPEng[i]);
        }
        if (scatpMomR.Mag() == 0.) 
        {
            for (unsigned int i = 0; i < OffMEng.GetSize(); i++)
            {
                if (OffMEng[i] < protonEnergy*0.5) continue; // Only look at protons with at least 50% of the beam energy
                scatpMomR = scatpMomT;
                scatp4MomR = scatp4MomT;
                break;
            }
            if (scatpMomR.Mag() == 0.)
            {
                for (unsigned int i = 0; i < B0Eng.GetSize(); i++)
                {
                    if (B0Eng[i] > 2) // Only look at tracks with at least 2 GeV energy
                    {
                        scatpMomR = scatpMomT;
                        scatp4MomR = scatp4MomT;
                    }
                }
            }
        } 

        if (scatpMomR.Mag() == 0.) // If no proton found, skip event
        { 
            cutEvents[1]++;
            continue; 
        }

        // Fill reconstructed track momentum vectors and find true PID
        for (unsigned int i = 0; i < trackEng.GetSize(); i++)
        {
            recoTrackMom[i] = TVector3(trackMomX[i],trackMomY[i],trackMomZ[i]);
            recoTrack4Mom[i].SetPxPyPzE(trackMomX[i],trackMomY[i],trackMomZ[i], TMath::Sqrt(trackMomX[i]*trackMomX[i] + trackMomY[i]*trackMomY[i] + trackMomZ[i]*trackMomZ[i] + muMass*muMass));
        
            /*
            if (TMath::Abs(recoTrackMom[i].Mag() - scatEMomT.Mag()) < 0.5)
            {
                recoTrackTruePID[i] = 11;     
            }
            else if (TMath::Abs(recoTrackMom[i].Mag() - muPlusMomT.Mag())< 0.5 || TMath::Abs(recoTrackMom[i].Mag() - muMinusMomT.Mag()) < 0.5)
            {
                recoTrackTruePID[i] = 13;     
            }
            else
            {
                recoTrackTruePID[i] = -1;
            }
            */

            if(partGenStat[simuAssoc[i]] == 1) // Select stable thrown particles
            {
                recoTrackTruePID[i] = TMath::Abs(partPdg[simuAssoc[i]]);
            }
            else
            {
                int motherID = ParentIndex(simuAssoc[i], partGenStat, partParI, partPdg, partParb, partPare);
                if (motherID == -1)
                {
                    recoTrackTruePID[i] = -1;
                }
                else recoTrackTruePID[i] = TMath::Abs(partPdg[motherID]);
            }
        }
        if ((recoTrackTruePID[0] == -1 || recoTrackTruePID[1] == -1 || recoTrackTruePID[2] == -1)) //|| (recoTrackTruePID[0] + recoTrackTruePID[1] + recoTrackTruePID[2] != 11 + 13 + 13))
        {
            incorrectParentIDrate++;
            cutEvents[2]++;
            continue;
        }

        // Search for the electron reconstructed tracks using the electron finder
        for (unsigned int i = 0; i < trackEng.GetSize(); i++)
        {
            if (trackCharge[i] == 1.0)
            {
                recoTrackIndex[1] = i;
                continue; // Only look at negative charged tracks
            }
            //recoTrackMom = TVector3(trackMomX[i],trackMomY[i],trackMomZ[i]);
            int simuID = simuAssoc[i];
            electronFinderReturn =  IsElectron(recoTrackMom[i], simuID, EcalBarrelEng, EcalEndcapPEng, EcalEndcapNEng, simuAssocEcalBarrel, simuAssocEcalEndcapP, simuAssocEcalEndcapN);
            

            if (electronFinderReturn == true && electronFound == true) // If more than one electron found, skip event
            {
                electronFound = false;
                break;
            }
            else if (electronFinderReturn == true && electronFound == false) // If electron found and no electron has been found yet
            {
                recoTrackIndex[0] = i;
                recoTrackIndex[2] = 3 - i - recoTrackIndex[1];
                electronFound = true;
            }
        }

        if (electronFound == true)
        {
            if (recoTrackTruePID[recoTrackIndex[0]] != 11) electronFinderIncorrectRate++;
            if (recoTrackTruePID[recoTrackIndex[0]] == 11) electronFinderCorrectRate++;
        } 
        else
        {
            
            cutEvents[3]++;
            continue;
        }

        for (unsigned int i = 0; i < trackEng.GetSize(); i++)
        {
            if (i == recoTrackIndex[0])
            {
                continue; // Ignore electron track
            }
            //recoTrackMom = TVector3(trackMomX[i],trackMomY[i],trackMomZ[i]);
            int simuID = simuAssoc[i];
            muonFinderReturn =  IsMuon(recoTrackMom[i], simuID, EcalBarrelEng, EcalEndcapPEng, EcalEndcapNEng, HcalBarrelEng, HcalEndcapPEng, HcalEndcapNEng, 
                                        simuAssocEcalBarrel, simuAssocEcalEndcapP, simuAssocEcalEndcapN, simuAssocHcalBarrel, simuAssocHcalEndcapP, simuAssocHcalEndcapN);
                

            if (muonFinderReturn == true && trackCharge[i] == 1.0)
            {
                recoTrackIndex[1] = i;
            }
            else if (muonFinderReturn == true && trackCharge[i] == -1.0) 
            {
                recoTrackIndex[2] = i;
            }
        }

        if (recoTrackIndex[1] != -1 && recoTrackIndex[2] != -1)
        {
            if (recoTrackTruePID[recoTrackIndex[1]] != 13) muonFinderIncorrectRate++;
            if (recoTrackTruePID[recoTrackIndex[1]] == 13) muonFinderCorrectRate++;
            if (recoTrackTruePID[recoTrackIndex[2]] != 13) muonFinderIncorrectRate++;
            if (recoTrackTruePID[recoTrackIndex[2]] == 13) muonFinderCorrectRate++;
        }
        else
        {
            
            cutEvents[4]++;
            continue;
        }

        if ((recoTrackIndex[0] != -1 && recoTrackIndex[1] != -1 && recoTrackIndex[2] != -1)) // If all tracks identified, fill momentum vectors
        {
            scatEMomR = recoTrackMom[recoTrackIndex[0]];
            recoTrack4Mom[recoTrackIndex[0]].SetPxPyPzE(trackMomX[recoTrackIndex[0]],trackMomY[recoTrackIndex[0]],trackMomZ[recoTrackIndex[0]], TMath::Sqrt(trackMomX[recoTrackIndex[0]]*trackMomX[recoTrackIndex[0]] + trackMomY[recoTrackIndex[0]]*trackMomY[recoTrackIndex[0]] + trackMomZ[recoTrackIndex[0]]*trackMomZ[recoTrackIndex[0]] + eMass*eMass));
            scatE4MomR = recoTrack4Mom[recoTrackIndex[0]];

            muPlusMomR = recoTrackMom[recoTrackIndex[1]];
            muPlus4MomR = recoTrack4Mom[recoTrackIndex[1]];
            muMinusMomR = recoTrackMom[recoTrackIndex[2]];
            muMinus4MomR = recoTrack4Mom[recoTrackIndex[2]];

            if (TMath::Abs(recoTrackTruePID[recoTrackIndex[0]]) != 11)
            { 
                incorrectElectronIDrate++;
            }

            if (TMath::Abs(recoTrackTruePID[recoTrackIndex[1]]) != 13)
            {
                incorrectMuonIDrate++;
            }

            if (TMath::Abs(recoTrackTruePID[recoTrackIndex[2]]) != 13)
            {
                incorrectMuonIDrate++;
            } 

        }
        else
        {
            cutEvents[4]++;
            continue;
        }

        JPsiMomR = muPlusMomR + muMinusMomR; 
        JPsi4MomR = muPlus4MomR + muMinus4MomR;


        if (JPsi4MomR.M() < 2.5 || JPsi4MomR.M() > 4.0) // If invariant mass of muon pair is outside J/Psi mass window, skip event
        { 
            cutEvents[5]++;
            continue;
        }

        ROOT::Math::Boost boostToJPsiRest(-JPsi4MomR.X()/JPsi4MomR.E(), -JPsi4MomR.Y()/JPsi4MomR.E(), -JPsi4MomR.Z()/JPsi4MomR.E());

        muPlus4MomR_Boosted = boostToJPsiRest(muPlus4MomR);
        muMinus4MomR_Boosted = boostToJPsiRest(muMinus4MomR);

        // Convert the boosted 4-vectors' spatial components to TVector3 and compute the angle
        TVector3 muPlusBoostVec(muPlus4MomR_Boosted.X(), muPlus4MomR_Boosted.Y(), muPlus4MomR_Boosted.Z());
        TVector3 muMinusBoostVec(muMinus4MomR_Boosted.X(), muMinus4MomR_Boosted.Y(), muMinus4MomR_Boosted.Z());
        double muonAngle = muPlusBoostVec.Angle(muMinusBoostVec);
        muonScatterAngle->Fill(muonAngle);

        if (TMath::Abs(TMath::Pi() - muonAngle) > 0.5) // If the angle between the two muons in the J/Psi rest frame is too large, skip event
        { 
            cutEvents[5]++;
            continue;
        }

        double missingMass = (beamE4Mom + beamp4Mom - scatE4MomR - scatp4MomR - muPlus4MomR - muMinus4MomR).M();
        double mm2 = (beamE4Mom + beamp4Mom - scatE4MomR - scatp4MomR - muPlus4MomR - muMinus4MomR).M2();
        matchedMissingMassEpToJPsiX->Fill(missingMass);
        matchedMissingMass2EpToJPsiX->Fill(mm2);
        if (missingMass > 2.0 || mm2 < -1) // If missing mass is too large, skip event
        { 
            cutEvents[5]++;
            continue; 
        }
        
        
        eventsPassed++;

        // Fill plots

        matchedElectronEta->Fill(scatEMomR.PseudoRapidity());
        electronEff->Fill(scatEMomR.PseudoRapidity());
        matchedElectronMomHist->Fill(scatEMomR.Mag());
        electronMomEff->Fill(scatEMomR.Mag());
        matchedElectronEtaMom->Fill(scatEMomR.PseudoRapidity(),scatEMomR.Mag());
        float deltaEta_electron = scatEMomT.PseudoRapidity() - scatEMomR.PseudoRapidity();
        float deltaPhi_electron = TVector2::Phi_mpi_pi(scatEMomT.Phi() - scatEMomR.Phi());
        float deltaR_electron = TMath::Sqrt(deltaEta_electron*deltaEta_electron + deltaPhi_electron*deltaPhi_electron);
        matchedElectronTrackDeltaR->Fill(deltaR_electron);
      

        matchedProtonEta->Fill(scatpMomR.PseudoRapidity());
        protonEff->Fill(scatpMomR.PseudoRapidity());
        matchedProtonMomHist->Fill(scatpMomR.Mag());
        protonMomEff->Fill(scatpMomR.Mag());
        matchedProtonEtaMom->Fill(scatpMomR.PseudoRapidity(),scatpMomR.Mag());
        matchedProtonEMinusPzHist->Fill(scatp4MomR.E() - scatpMomR.Z());
        matchedProtonPtHist->Fill(scatpMomR.Perp());
        float deltaEta_proton = scatpMomT.PseudoRapidity() - scatpMomR.PseudoRapidity();
        float deltaPhi_proton = TVector2::Phi_mpi_pi(scatpMomT.Phi() - scatpMomR.Phi());
        float deltaR_proton = TMath::Sqrt(deltaEta_proton*deltaEta_proton + deltaPhi_proton*deltaPhi_proton);
        matchedProtonTrackDeltaR->Fill(deltaR_proton);
        
        matchedMuonEta->Fill(muPlusMomR.PseudoRapidity());
        muonEff->Fill(muPlusMomR.PseudoRapidity());
        matchedMuonMomHist->Fill(muPlusMomR.Mag());
        muonMomEff->Fill(muPlusMomR.Mag());
        matchedMuonEtaMom->Fill(muPlusMomR.PseudoRapidity(),muPlusMomR.Mag());
        float deltaEta_muPlus = muPlusMomT.PseudoRapidity() - muPlusMomR.PseudoRapidity();
        float deltaPhi_muPlus = TVector2::Phi_mpi_pi(muPlusMomT.Phi() - muPlusMomR.Phi());
        float deltaR_muPlus = TMath::Sqrt(deltaEta_muPlus*deltaEta_muPlus + deltaPhi_muPlus*deltaPhi_muPlus);
        matchedMuonTrackDeltaR->Fill(deltaR_muPlus);
        
        matchedElectronProtonEMinusPzHist->Fill(scatE4MomR.E() + scatp4MomR.E() - (scatE4MomR.Pz() + scatp4MomR.Pz()));
        matchedElectronProtonPtHist->Fill((scatEMomR + scatpMomR).Perp());
        
        matchedMuonEta->Fill(muMinusMomR.PseudoRapidity());
        muonEff->Fill(muMinusMomR.PseudoRapidity());
        matchedMuonMomHist->Fill(muMinusMomR.Mag());
        muonMomEff->Fill(muMinusMomR.Mag());
        matchedMuonEtaMom->Fill(muMinusMomR.PseudoRapidity(),muMinusMomR.Mag()); 
        float deltaEta_muMinus = muMinusMomT.PseudoRapidity() - muMinusMomR.PseudoRapidity();
        float deltaPhi_muMinus = TVector2::Phi_mpi_pi(muMinusMomT.Phi() - muMinusMomR.Phi());
        float deltaR_muMinus = TMath::Sqrt(deltaEta_muMinus*deltaEta_muMinus + deltaPhi_muMinus*deltaPhi_muMinus);
        matchedMuonTrackDeltaR->Fill(deltaR_muMinus);

        matchedMuMuInvMass->Fill(JPsi4MomR.M()); 

        matchedJPsiEta->Fill(JPsiMomR.Eta());
        JPsiEff->Fill(JPsiMomR.Eta());
        matchedJPsiMomHist->Fill(JPsiMomR.Mag());
        JPsiMomEff->Fill(JPsiMomR.Mag());
        matchedJPsiEtaMom->Fill(JPsiMomR.Eta(), JPsiMomR.Mag());  
        float deltaEta_JPsi = JPsiMomT.Eta() - JPsiMomR.Eta();
        float deltaPhi_JPsi = TVector2::Phi_mpi_pi(JPsiMomT.Phi() - JPsiMomR.Phi());
        float deltaR_JPsi = TMath::Sqrt(deltaEta_JPsi*deltaEta_JPsi + deltaPhi_JPsi*deltaPhi_JPsi);
        matchedJPsiTrackDeltaR->Fill(deltaR_JPsi);       

        // Calculate kinematic variable

        Q2_truth = -(beamE4Mom - scatE4MomT).mag2();
        t_truth = -1*((scatp4MomT - beamp4Mom).mag2());
        y_truth =(beamp4Mom.Dot(beamE4Mom - scatE4MomT))/(beamp4Mom.Dot(beamE4Mom));
        x_truth = Q2_truth/(4*beamE4Mom.E()*beamp4Mom.E()*y_truth);

        double delta_h = (JPsi4MomR.E() + scatp4MomR.E()) - (JPsi4MomR.Pz() + scatp4MomR.Pz());
        double pt2_h = (pow((JPsi4MomR.Px()+scatp4MomR.Px()),2))+(pow((JPsi4MomR.Py()+scatp4MomR.Py()),2));
        double alpha_e = tan((scatE4MomR.Theta()/2));
        double alpha_h = delta_h/(sqrt(pt2_h));

        ROOT::Math::PxPyPzEVector vec_t_XPT_sum = (scatE4MomR + JPsi4MomR);
        TVector3 vec_t_XPT3 = TVector3(vec_t_XPT_sum.Px(), vec_t_XPT_sum.Py(), vec_t_XPT_sum.Pz());
        vec_t_XPT3.RotateY(crossingAngle);
        ROOT::Math::PxPyPzEVector vec_t_XPT = ROOT::Math::PxPyPzEVector(vec_t_XPT3.X(), vec_t_XPT3.Y(), vec_t_XPT3.Z(), vec_t_XPT_sum.E());

        ROOT::Math::PxPyPzEVector vec_PMiss_Rec = (beamE4Mom + beamp4Mom) - (scatE4MomR + JPsi4MomR);
        ROOT::Math::PxPyPzEVector scatp4MomR_corr = ROOT::Math::PxPyPzEVector(vec_PMiss_Rec.P()*sin(scatp4MomT.Theta())*cos(scatp4MomT.Phi()), vec_PMiss_Rec.P()*sin(scatp4MomT.Theta())*sin(scatp4MomT.Phi()), vec_PMiss_Rec.P()*cos(scatp4MomT.Theta()), sqrt(pow(vec_PMiss_Rec.P(),2)+(pow(pMass,2))));

        y_e =(beamp4Mom.Dot(beamE4Mom - scatE4MomR))/(beamE4Mom.Dot(beamE4Mom));
        Q2_e = -(beamE4Mom - scatE4MomR).mag2();
        x_e = Q2_e/(4*beamE4Mom.E()*beamp4Mom.E()*y_e);

        y_JB = delta_h/(2*beamE4Mom.E());
        Q2_JB = pt2_h/(1-y_JB);
        x_JB = Q2_JB/(4*beamE4Mom.E()*beamp4Mom.E()*y_JB);

        y_DA = (alpha_h)/(alpha_e + alpha_h);
        Q2_DA = (4*beamE4Mom.E()*beamE4Mom.E())/(alpha_e*(alpha_e +alpha_h));
        x_DA = Q2_DA/(4*beamE4Mom.E()*beamp4Mom.E()*y_DA);
        
        y_sigma = delta_h/(delta_h + (scatE4MomR.E()*(1-cos(scatE4MomR.Theta()))));
        Q2_sigma = (pow(scatE4MomR.E(),2)*pow(sin(scatE4MomR.Theta()),2))/(1-y_sigma);
        x_sigma = Q2_sigma/(4*beamE4Mom.E()*beamp4Mom.E()*y_sigma);

        t_BABE = -1*((beamp4Mom - scatp4MomR).mag2());
        t_eX = -1*(((beamE4Mom - scatE4MomR) - JPsi4MomR).mag2());
        t_eXPT = (vec_t_XPT.Perp2()); // Rotate vetors prior to getting perpendicular component
        t_eXBABE = -1*((beamp4Mom - scatp4MomR_corr).mag2());

        // Fill kinematic histograms

        trueQ2->Fill(Q2_truth);
        reconQ2_e->Fill(Q2_e);
        reconQ2_JB->Fill(Q2_JB);
        reconQ2_DA->Fill(Q2_DA);
        reconQ2_sigma->Fill(Q2_sigma);
        deltaQ2_e->Fill(100*(Q2_e - Q2_truth)/Q2_truth);
        deltaQ2_JB->Fill(100*(Q2_JB - Q2_truth)/Q2_truth);
        deltaQ2_DA->Fill(100*(Q2_DA - Q2_truth)/Q2_truth);
        deltaQ2_sigma->Fill(100*(Q2_sigma - Q2_truth)/Q2_truth);

        reconQ2_DA_vs_trueQ2->Fill(Q2_truth, Q2_DA);

        truet->Fill(t_truth);
        recont_eXBABE->Fill(t_eXBABE);
        recont_eXPT->Fill(t_eXPT);
        recont_eX->Fill(t_eX);
        recont_BABE->Fill(t_BABE);
        deltat_eXBABE->Fill(100*(t_eXBABE - t_truth)/t_truth);
        deltat_eXPT->Fill(100*(t_eXPT - t_truth)/t_truth);
        deltat_eX->Fill(100*(t_eX - t_truth)/t_truth);
        deltat_BABE->Fill(100*(t_BABE - t_truth)/t_truth);

        recont_eXBABE_vs_truet->Fill(t_truth, t_eXBABE);
        deltat_eXBABE_vs_truet->Fill(t_truth, 100*(t_eXBABE - t_truth)/t_truth);

        truey->Fill(y_truth);
        recony_e->Fill(y_e);
        recony_JB->Fill(y_JB);
        recony_DA->Fill(y_DA);
        recony_sigma->Fill(y_sigma);
        deltay_e->Fill(100*(y_e - y_truth)/y_truth);
        deltay_JB->Fill(100*(y_JB - y_truth)/y_truth);
        deltay_DA->Fill(100*(y_DA - y_truth)/y_truth);
        deltay_sigma->Fill(100*(y_sigma - y_truth)/y_truth);

        recony_DA_vs_truey->Fill(y_truth, y_DA);

        truex->Fill(x_truth);
        reconx_e->Fill(x_e);
        reconx_JB->Fill(x_JB);
        reconx_DA->Fill(x_DA);
        reconx_sigma->Fill(x_sigma);
        deltax_e->Fill(100*(x_e - x_truth)/x_truth);
        deltax_JB->Fill(100*(x_JB - x_truth)/x_truth);
        deltax_DA->Fill(100*(x_DA - x_truth)/x_truth);
        deltax_sigma->Fill(100*(x_sigma - x_truth)/x_truth);


        if (1 <= Q2_truth && Q2_truth <= 50) // Check if Q2 is in the range of interest
        {
          if (0.0016 <= x_truth && x_truth < 0.0025) // Check if xbjk is in the range of interest
          {
            truet_XbjkA->Fill(t_truth);
          }
          else if (0.016 <= x_truth && x_truth < 0.025)
          {
            truet_XbjkB->Fill(t_truth);
          }
          else if (0.16 <= x_truth && x_truth < 0.25)
          {
            truet_XbjkC->Fill(t_truth);
          }
        }

        if (1 <= Q2_DA && Q2_DA <= 50) // Check if Q2 is in the range of interest
        {
          if (0.0016 <= x_DA && x_DA < 0.0025) // Check if xbjk is in the range of interest
          {
            recont_XbjkA->Fill(t_eXBABE);
          }
          else if (0.016 <= x_DA && x_DA < 0.025)
          {
            recont_XbjkB->Fill(t_eXBABE);
          }
          else if (0.16 <= x_DA && x_DA < 0.25)
          {
            recont_XbjkC->Fill(t_eXBABE);
          }
        }
    } 

    std::cout << "Event Processing Complete" << std::endl;
    std::cout << "Total Events Processed: " << eventID << std::endl;
    std::cout << "Cut Flow: " << std::endl;
    std::cout << "  - Events with 3 tracks and sum charge of -1: " << eventID - cutEvents[0] << std::endl;
    std::cout << "  - Events with identified proton: " << eventID - cutEvents[0] - cutEvents[1] << std::endl;
    std::cout << "  - Events with correct parents: " << eventID - cutEvents[0] - cutEvents[1] - cutEvents[2] << std::endl;
    std::cout << "  - Events with identified electron: " << eventID - cutEvents[0] - cutEvents[1] - cutEvents[2] - cutEvents[3] << std::endl;
    std::cout << "  - Events with identified muon pair: " << eventID - cutEvents[0] - cutEvents[1] - cutEvents[2] - cutEvents[3] - cutEvents[4] << std::endl;
    std::cout << "  - Events with J/Psi candidate from muon pair: " << eventID - cutEvents[0] - cutEvents[1] - cutEvents[2] - cutEvents[3] - cutEvents[4] - cutEvents[5] << std::endl;
    std::cout << "  - Total Events Passing All Cuts: " << eventsPassed << std::endl;
    std::cout << "Parent ID Incorrect Rate: " << incorrectParentIDrate << std::endl;
    std::cout << "Electron Finder Correct Rate: " << electronFinderCorrectRate << ", Inconclusive Rate: " << electronFinderInconclusiveRate << ", Fail Rate: " << electronFinderIncorrectRate << std::endl;
    std::cout << "Muon Finder Correct Rate: " << muonFinderCorrectRate << ", Inconclusive Rate: " << muonFinderInconclusiveRate << ", Fail Rate: " << muonFinderIncorrectRate << std::endl;
    std::cout << "Incorrectly identified electron tracks: " << incorrectElectronIDrate << std::endl;
    std::cout << "Incorrectly identified muon tracks: " << incorrectMuonIDrate << std::endl;
    
    electronEff->Divide(electronEta);
    protonEff->Divide(protonEta);
    muonEff->Divide(muonEta);
    JPsiEff->Divide(JPsiEta);

    electronMomEff->Divide(electronMomHist);
    protonMomEff->Divide(protonMomHist);
    muonMomEff->Divide(muonMomHist);
    JPsiMomEff->Divide(JPsiMomHist);

    // Write histograms to file
    ofile->cd();
    ofile->mkdir("electrons");
    ofile->cd("electrons");
    electronEta->Write();
    matchedElectronEta->Write();
    electronEff->Write();
    electronMomHist->Write();
    matchedElectronMomHist->Write();
    electronMomEff->Write();
    electronEtaMom->Write();
    matchedElectronEtaMom->Write();
    matchedElectronTrackDeltaR->Write();
    ofile->cd("..");
    ofile->mkdir("protons");
    ofile->cd("protons");
    protonEta->Write();
    matchedProtonEta->Write();
    protonEff->Write();
    protonMomHist->Write();
    matchedProtonMomHist->Write();
    protonMomEff->Write();
    protonEtaMom->Write();
    protonEMinusPzHist->Write();
    matchedProtonEMinusPzHist->Write();
    protonPtHist->Write();
    matchedProtonPtHist->Write();
    matchedProtonEtaMom->Write();
    matchedProtonTrackDeltaR->Write();
    ofile->cd("..");
    ofile->mkdir("electronProton");
    ofile->cd("electronProton");
    electronProtonEMinusPzHist->Write();
    matchedElectronProtonEMinusPzHist->Write();
    electronProtonPtHist->Write();
    matchedElectronProtonPtHist->Write();
    ofile->cd("..");
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
    ofile->cd("..");
    ofile->mkdir("invariantMass");
    ofile->cd("invariantMass");
    muonScatterAngle->Write();
    invMassMuMu->Write();
    matchedMuMuInvMass->Write();
    matchedMissingMassEpToJPsiX->Write();
    matchedMissingMass2EpToJPsiX->Write();
    ofile->cd("..");
    ofile->mkdir("kinematics");
    ofile->cd("kinematics");
    trueQ2->Write();
    reconQ2_e->Write();
    reconQ2_JB->Write();
    reconQ2_DA->Write();
    reconQ2_sigma->Write();
    reconQ2_DA_vs_trueQ2->Write();
    truet->Write();
    recont_eXBABE->Write();
    recont_eXPT->Write();
    recont_eX->Write();
    recont_BABE->Write();
    recont_eXBABE_vs_truet->Write();
    truey->Write();
    recony_e->Write();
    recony_JB->Write();
    recony_DA->Write();
    recony_sigma->Write();
    recony_DA_vs_truey->Write();
    truex->Write();
    reconx_e->Write();
    reconx_JB->Write();
    reconx_DA->Write();
    reconx_sigma->Write();
    ofile->cd("..");
    ofile->mkdir("kinematicDifferences");
    ofile->cd("kinematicDifferences");
    deltaQ2_e->Write();
    deltaQ2_JB->Write();
    deltaQ2_DA->Write();
    deltaQ2_sigma->Write();
    deltat_eXBABE->Write();
    deltat_eXPT->Write();
    deltat_eX->Write();
    deltat_BABE->Write();
    deltat_eXBABE_vs_truet->Write();
    deltay_e->Write();
    deltay_JB->Write();
    deltay_DA->Write();
    deltay_sigma->Write();
    deltax_e->Write();
    deltax_JB->Write();
    deltax_DA->Write();
    deltax_sigma->Write();
    ofile->cd("..");
    ofile->mkdir("tDistributionsByXbjk");
    ofile->cd("tDistributionsByXbjk");
    truet_XbjkA->Write();
    truet_XbjkB->Write();
    truet_XbjkC->Write();
    recont_XbjkA->Write();
    recont_XbjkB->Write();
    recont_XbjkC->Write();
    ofile->cd("..");


    ofile->Close(); // Close output file

    std::cout << "Histograms written to file: " << outfilename << std::endl;
    std::cout << "Track analysis completed." << std::endl;

}