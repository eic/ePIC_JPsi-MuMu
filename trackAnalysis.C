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

void trackAnalysis()
  {

    TString infile="eicReconOutput/EICreconOut_JPsiMuMu_10ifb_18x275ep_Pruned.root";
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
    std::string outfilename = "outputs/trackAnalysisOutput_" + std::to_string(lumi_int) + "ifb_" + beam_config + ".root";

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

    // Get Forward Detector Information
    TTreeReaderArray<float> RPEng(tree_reader, "ForwardRomanPotRecParticles.energy");
    TTreeReaderArray<float> RPMomX(tree_reader, "ForwardRomanPotRecParticles.momentum.x");
    TTreeReaderArray<float> RPMomY(tree_reader, "ForwardRomanPotRecParticles.momentum.y");
    TTreeReaderArray<float> RPMomZ(tree_reader, "ForwardRomanPotRecParticles.momentum.z");
    TTreeReaderArray<float> OffMEng(tree_reader, "ForwardOffMRecParticles.energy");

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
    TH1D *matchedElectronMomHist = new TH1D("matchedElectronMomHist","Momentum of Thrown Electrons That Have Matching Track;Momentum (GeV/c)",100,-100.,100.);
    TH1D *electronMomEff = new TH1D("electronMomEff","Efficency;Momentum (GeV/c)",100,-100.,100.);

    TH1D *protonMomHist = new TH1D("protonMomHist","Momentum of Thrown Protons;Momentum (GeV/c)",160,-20.,300);
    TH1D *matchedProtonMomHist = new TH1D("matchedProtonMomHist","Momentum of Thrown Protons That Have Matching Track;Momentum (GeV/c)",160,-20.,300);
    TH1D *protonMomEff = new TH1D("protonMomEff","Efficency;Momentum (GeV/c)",160,-20.,300);

    TH1D *muonMomHist = new TH1D("muonMomHist","Momentum of Thrown Muons;Momentum (GeV/c)",160,-20.,300);
    TH1D *matchedMuonMomHist = new TH1D("matchedMuonMomHist","Momentum of Thrown Muons That Have Matching Track;Momentum (GeV/c)",160,-20.,300);
    TH1D *muonMomEff = new TH1D("muonMomEff","Efficency;Momentum (GeV/c)",160,-20.,300);

    TH1D *JPsiMomHist = new TH1D("JPsiMomHist","Momentum of Thrown J/PSi;Momentum (GeV/c)",160,-20.,300);
    TH1D *matchedJPsiMomHist = new TH1D("matchedJPsiMomHist","Momentum of Thrown J/PSi That Have Matching Track;Momentum (GeV/c)",160,-20.,300);
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

    // Define Invariant Mass Histograms

    TH1D *invMassMuMu = new TH1D("invMassMuMu","Invariant Mass of Muon Pairs",100,0.,5.); // Invariant mass histogram
    TH1D *matchedMuMuInvMass = new TH1D("matchedMuMuInvMass","Invariant Mass of Muon Pairs That Have Matching Track",100,0.,5.); // Invariant mass histogram

    // Define Kinematics Histograms
    TH1D* trueQ2 = new TH1D("trueQ^{2}","True Q^{2} Distribution; Q^{2} (GeV/c^{2})",100,0.,50.);
    TH1D* reconQ2_DA = new TH1D("reconQ^{2}_DA","Reconstructed Q^{2} Distribution using the DA Method; Q^{2} (GeV/c^{2})",100,0.,50.);
    TH1D* reconQ2_JB = new TH1D("reconQ^{2}_JB","Reconstructed Q^{2} Distribution using the JB Method; Q^{2} (GeV/c^{2})",100,0.,50.);
    TH1D* reconQ2_e = new TH1D("reconQ^{2}_e","Reconstructed Q^{2} Distribution using the electron Method; Q^{2} (GeV/c^{2})",100,0.,50.);
    TH1D* reconQ2_sigma = new TH1D("reconQ^{2}_sigma","Reconstructed Q^{2} Distribution using the sigma Method; Q^{2} (GeV/c^{2})",100,0.,50.);

    TH1D* deltaQ2_DA = new TH1D("deltaQ2_DA","Delta Q^{2} (Reconstructed - True) using the DA Method; (Q^{2} - Q_{MC}^{2})/Q_{MC}^{2}",100,-10.,10.);
    TH1D* deltaQ2_JB = new TH1D("deltaQ2_JB","Delta Q^{2} (Reconstructed - True) using the JB Method; (Q^{2} - Q_{MC}^{2})/Q_{MC}^{2}",100,-10.,10.);
    TH1D* deltaQ2_e = new TH1D("deltaQ2_e","Delta Q^{2} (Reconstructed - True) using the electron Method; (Q^{2} - Q_{MC}^{2})/Q_{MC}^{2}",100,-10.,10.);
    TH1D* deltaQ2_sigma = new TH1D("deltaQ2_sigma","Delta Q^{2} (Reconstructed - True) using the sigma Method; (Q^{2} - Q_{MC}^{2})/Q_{MC}^{2}",100,-10.,10.);

    TH2D* reconQ2_DA_vs_trueQ2 = new TH2D("reconQ2_DA_vs_trueQ2","Reconstructed Q^{2} using the DA Method vs True Q^{2}; Q_{MC}^{2} (GeV/c^{2}); Q_{DA}^{2} (GeV/c^{2})",100,0.,50.,100,0.,50.);

    TH1D* truet = new TH1D("truet","True t Distribution; t (GeV/c)",100,0.,2.);
    TH1D* recont_eXBABE = new TH1D("recont_eXBABE","Reconstructed t Distribution using the eXBABE Method; t (GeV/c)",100,0.,2.);
    TH1D* recont_eXPT = new TH1D("recont_eXPT","Reconstructed t Distribution using the eXPT Method; t (GeV/c)",100,0.,2.);
    TH1D* recont_eX = new TH1D("recont_eX","Reconstructed t Distribution using the eX Method; t (GeV/c)",100,0.,2.);
    TH1D* recont_BABE = new TH1D("recont_BABE","Reconstructed t Distribution using the BABE Method; t (GeV/c)",100,0.,2.);

    TH1D* deltat_eXBABE = new TH1D("deltat_eXBABE","Delta t (Reconstructed - True) using the eXBABE Method; (t - t_{MC})/t_{MC}",100,-0.5,0.5);
    TH1D* deltat_eXPT = new TH1D("deltat_eXPT","Delta t (Reconstructed - True) using the eXPT Method; (t - t_{MC})/t_{MC}",100,-0.5,0.5);
    TH1D* deltat_eX = new TH1D("deltat_eX","Delta t (Reconstructed - True) using the eX Method; (t - t_{MC})/t_{MC}",100,-0.5,0.5);
    TH1D* deltat_BABE = new TH1D("deltat_BABE","Delta t (Reconstructed - True) using the BABE Method; (t - t_{MC})/t_{MC}",100,-0.5,0.5);

    TH2D* recont_eXBABE_vs_truet = new TH2D("recont_eXBABE_vs_truet","Reconstructed t using the eXBABE Method vs True t; t_{MC} (GeV/c); t_{eXBABE} (GeV/c)",100,0.,2.,100,0.,2.);

    TH1D* truex = new TH1D("truex","True x Distribution",100,0.,0.25);
    TH1D* reconx_DA = new TH1D("reconx_DA","Reconstructed x Distribution using the DA method; x",100,0.,0.25);
    TH1D* reconx_JB = new TH1D("reconx_JB","Reconstructed x Distribution using the JB method; x",100,0.,0.25);
    TH1D* reconx_e = new TH1D("reconx_e","Reconstructed x Distribution using the electron method; x",100,0.,0.25);
    TH1D* reconx_sigma = new TH1D("reconx_sigma","Reconstructed x Distribution using the sigma method; x",100,0.,0.25);

    TH1D* deltax_DA = new TH1D("deltax_DA","Delta x (Reconstructed - True) using the DA Method; (x - x_{MC})/x_{MC}",100,-1.,1.);
    TH1D* deltax_JB = new TH1D("deltax_JB","Delta x (Reconstructed - True) using the JB Method; (x - x_{MC})/x_{MC}",100,-1.,1.);
    TH1D* deltax_e = new TH1D("deltax_e","Delta x (Reconstructed - True) using the electron Method; (x - x_{MC})/x_{MC}",100,-1.,1.);
    TH1D* deltax_sigma = new TH1D("deltax_sigma","Delta x (Reconstructed - True) using the sigma Method; (x - x_{MC})/x_{MC}",100,-1.,1.);

    TH1D* truet_XbjkA = new TH1D("truet_XbjkA","True t distribution with bjorken x binning",100,0.,2.);
    TH1D* recont_XbjkA = new TH1D("recont_XbjkA","Reconstructed t distribution with bjorken x binning",100,0.,2.);
    TH1D* truet_XbjkB = new TH1D("truet_XbjkB","True t distribution with bjorken x binning",100,0.,2.);
    TH1D* recont_XbjkB = new TH1D("recont_XbjkB","Reconstructed t distribution with bjorken x binning",100,0.,2.);
    TH1D* truet_XbjkC = new TH1D("truet_XbjkC","True t distribution with bjorken x binning",100,0.,2.);
    TH1D* recont_XbjkC = new TH1D("recont_XbjkC","Reconstructed t distribution with bjorken x binning",100,0.,2.);

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

    double Q2_truth, Q2_DA, Q2_JB, Q2_e, Q2_sigma;
    double t_truth, t_eXBABE, t_eXPT, t_eX, t_BABE;
    double x_truth, x_DA, x_JB, x_e, x_sigma;
    double y_truth, y_DA, y_JB, y_e, y_sigma;
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

      // Reset kinematic variables
      Q2_truth = 0.;
      Q2_DA = 0.;
      Q2_JB = 0.;
      Q2_e = 0.;
      Q2_sigma = 0.;
      t_truth = 0.;
      t_eXBABE = 0.;
      t_eXPT = 0.;
      t_eX = 0.;
      t_BABE = 0.;
      x_truth = 0.;
      x_DA = 0.;
      x_JB = 0.;
      x_e = 0.;
      x_sigma = 0.;
      y_truth = 0.;
      y_DA = 0.;
      y_JB = 0.;
      y_e = 0.;
      y_sigma = 0.;

      for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop over thrown particles
        {
          if(partGenStat[i] == 1) // Select stable thrown particles
            {
              int pdg = TMath::Abs(partPdg[i]);
              partEng = sqrt(pow(partMomX[i],2) + pow(partMomY[i],2) + pow(partMomZ[i],2) + pow(partMass[i],2)); // Energy of all Monte Carlo particles

              if(pdg == 11) // Look at electrons
              {
                scatEMomT = TVector3(partMomX[i],partMomY[i],partMomZ[i]);
                scatE4MomT.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng); 
              
                // Loop over associations to find matching ReconstructedChargedParticle
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
                      scatEMomR = TVector3(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle
                      scatE4MomR.SetPxPyPzE(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]], trackEng[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle
                      break;

                    }
                  }
                }
                if (scatE == false)
                {
                  for (int iTtrack = 0; iTtrack < ttrackMomX.GetSize(); iTtrack++) // Look for far-backwards tracks
                  {
                    TVector3 trecMom(ttrackMomX[iTtrack],ttrackMomY[iTtrack],ttrackMomZ[iTtrack]);
                    ROOT::Math::PxPyPzEVector trec4Mom(ttrackMomX[iTtrack],ttrackMomY[iTtrack],ttrackMomZ[iTtrack], ttrackEng[iTtrack]); 

                    if ((ttrackPDG[iTtrack] == 0 || ttrackPDG[iTtrack] == 11) && trec4Mom.Eta() < -4.0 && ttrackCharge[iTtrack] == -1)
                    {
                      scatEMomR = trecMom;
                      scatE4MomR = trec4Mom;
                      scatE = true;
                      break; 
                    }
                  }              
                }    
              } 
              if(pdg == 13 && partCharge[i] == -1) // Look at muons
                {
                  muMinusMomT = TVector3(partMomX[i],partMomY[i],partMomZ[i]);
                  muMinus4MomT.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);

                  // Loop over associations to find matching ReconstructedChargedParticle
                  for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
                    {
                      if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
                        {
                          muMinusMomR = TVector3(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle
                          muMinus4MomR.SetPxPyPzE(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]], trackEng[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle
                          
                        }
                      else
                        {
                          for (int iTtrack = 0; iTtrack < ttrackMomX.GetSize(); iTtrack++) // Look for far-forwards tracks
                          {
                            TVector3 trecMom(ttrackMomX[iTtrack],ttrackMomY[iTtrack],ttrackMomZ[iTtrack]);
                            ROOT::Math::PxPyPzEVector trec4Mom(ttrackMomX[iTtrack],ttrackMomY[iTtrack],ttrackMomZ[iTtrack], ttrackEng[iTtrack]); 

                            if ((ttrackPDG[iTtrack] == 13) && ttrackCharge[iTtrack] == -1)
                            {
                              muMinusMomR = trecMom;
                              muMinus4MomR = trec4Mom;
                              break; 

                            }
                          }                   
                        }
                    }
                } 
              if(pdg == 13 && partCharge[i] == 1) // Look at muons
                {
                  muPlusMomT = TVector3(partMomX[i],partMomY[i],partMomZ[i]);
                  muPlus4MomT.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);

                  // Loop over associations to find matching ReconstructedChargedParticle
                  for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
                    {
                      if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
                        {
                          muPlusMomR = TVector3(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle
                          muPlus4MomR.SetPxPyPzE(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]], trackEng[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle
                          
                        }
                      else
                        {
                          for (int iTtrack = 0; iTtrack < ttrackMomX.GetSize(); iTtrack++) // Look for far-forwards tracks
                          {
                            TVector3 trecMom(ttrackMomX[iTtrack],ttrackMomY[iTtrack],ttrackMomZ[iTtrack]);
                            ROOT::Math::PxPyPzEVector trec4Mom(ttrackMomX[iTtrack],ttrackMomY[iTtrack],ttrackMomZ[iTtrack], ttrackEng[iTtrack]); 

                            if ((ttrackPDG[iTtrack] == 13) && ttrackCharge[iTtrack] == 1)
                            {
                              muPlusMomR = trecMom;
                              muPlus4MomR = trec4Mom;
                              break; 

                            }
                          }                   
                        }
                    }
                }
              if(pdg == 2212)
              {
                scatpMomT = TVector3(partMomX[i],partMomY[i],partMomZ[i]);
                scatpMomT.RotateY(0.025);
                scatp4MomT.SetPxPyPzE(scatpMomT.X(), scatpMomT.Y(), scatpMomT.Z() , partEng);

                // Loop over associations to find matching ReconstructedChargedParticle
                for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
                  {
                    if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
                      {
                        //std::cout << "Proton found with energy: " << partEng << " GeV" << std::endl;
                        scatpMomR = TVector3(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle
                        scatp4MomR.SetPxPyPzE(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]], trackEng[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle

                      }
                  }
                  if (scatpMomR.Mag() > 50) continue;
                    
                      for (int iTtrack = 0; iTtrack < ttrackMomX.GetSize(); iTtrack++) // Look for far-forwards tracks
                      {
                        TVector3 trecMom(ttrackMomX[iTtrack],ttrackMomY[iTtrack],ttrackMomZ[iTtrack]);
                        ROOT::Math::PxPyPzEVector trec4Mom(ttrackMomX[iTtrack],ttrackMomY[iTtrack],ttrackMomZ[iTtrack], ttrackEng[iTtrack]); 

                        if ((ttrackPDG[iTtrack] == 0 || ttrackPDG[iTtrack] == 2212) && trec4Mom.Eta() >4.0 && ttrackCharge[iTtrack] == 1 && ttrackEng[iTtrack] > 50)
                        {
                          //std::cout << "Found positive particle in far-forward with energy: " << ttrackEng[iTtrack] << " GeV" << std::endl;
                          scatpMomR = trecMom;
                          scatp4MomR = trec4Mom;
                          break;

                        }
                      }   

                  if (scatpMomR.Mag() > 50) continue;
                    for (int k = 0; k < B0Eng.GetSize(); k++)
                      {
                        if(B0Eng[k] > 50) // Look for hit in the B0
                        {
                          //std::cout << "Found particle in B0 with energy: " << B0Eng[k] << " GeV" << std::endl;
                          scatpMomR = scatpMomT; 
                          scatp4MomR = scatp4MomT; 
                        }       
                      }  

                  if (scatpMomR.Mag() > 50) continue; //std::cout << "Proton not found in truth-seeded far-forward tracks" << std::endl;              
                    for (int r = 0; r < RPEng.GetSize(); r++) // Look in the Roman Pots
                      {
                        if(RPEng[r] > 50) // Look for hit in the RP
                        {
                          //std::cout << "Found particle in RP with energy: " << RPEng[r] << " GeV" << std::endl;
                          scatpMomR = scatpMomT; 
                          scatp4MomR = scatp4MomT;
                          break;
                        }       
                      }
                      //std::cout << "Proton not found in B0" << std::endl;
                  if (scatpMomR.Mag() > 50) continue;
                  for (int m = 0; m < OffMEng.GetSize(); m++) // Look in the Off-Momentum Detectors
                    {
                      if(OffMEng[m] > 50) // Look for hit in the Off-Momentum Detectors
                      {
                        //std::cout << "Found particle in Off-Momentum Detectors with energy: " << OffMEng[m] << " GeV" << std::endl;
                        scatpMomR = scatpMomT; 
                        scatp4MomR = scatp4MomT; 
                      }     
                    }   
                    
              }
          }    
        }

      // Process the event

      if (scatEMomT.Mag() > 0 && scatpMomT.Mag() > 0 && muPlusMomT.Mag() > 0 && muMinusMomT.Mag() > 0) // Check if the scattered electron, proton, and muons have been found
      {

        electronEta->Fill(scatEMomT.PseudoRapidity());
        electronMomHist->Fill(scatEMomT.Mag());
        electronEtaMom->Fill(scatEMomT.PseudoRapidity(),scatEMomT.Mag());
      
        protonEta->Fill(scatpMomT.PseudoRapidity());
        protonMomHist->Fill(scatpMomT.Mag());
        protonEtaMom->Fill(scatpMomT.PseudoRapidity(),scatpMomT.Mag());
      
        muonEta->Fill(muPlusMomT.PseudoRapidity());
        muonMomHist->Fill(muPlusMomT.Mag());
        muonEtaMom->Fill(muPlusMomT.PseudoRapidity(),muPlusMomT.Mag());
        
        muonEta->Fill(muMinusMomT.PseudoRapidity());
        muonMomHist->Fill(muMinusMomT.Mag());
        muonEtaMom->Fill(muMinusMomT.PseudoRapidity(),muMinusMomT.Mag()); 

        JPsiMomT = muPlusMomT + muMinusMomT; 
        JPsi4MomT = muPlus4MomT + muMinus4MomT; 

        invMassMuMu->Fill(JPsi4MomT.M());

        JPsiEta->Fill(JPsiMomT.Eta());
        JPsiMomHist->Fill(JPsiMomT.Mag());
        JPsiEtaMom->Fill(JPsiMomT.Eta(), JPsiMomT.Mag());


        Q2_truth = -(beamE4Mom - scatE4MomT).mag2();
        t_truth = -1*((scatp4MomT - beamp4Mom).mag2());
        y_truth =(beamp4Mom.Dot(beamE4Mom - scatE4MomT))/(beamp4Mom.Dot(beamE4Mom));
        x_truth = Q2_truth/(4*beamE4Mom.E()*beamp4Mom.E()*y_truth);
        
        //std::cout << "Event ID: " << eventID << " Q2: " << Q2_truth << " t: " << t_truth << " xbjk: " << x_truth << std::endl;
        
        trueQ2->Fill(Q2_truth);
        truet->Fill(t_truth); 
        truex->Fill(x_truth);

        // Check if the scattered electron, proton, and muons have been found in the reconstructed tracks
        //std::cout << "Scattered Electron: " << scatEMomR.Mag() << " Scattered Proton: " << scatpMomR.Mag() << " Mu+ Momentum: " << muPlusMomR.Mag() << " Mu- Momentum: " << muMinusMomR.Mag() << std::endl;
    
        if (scatEMomR.Mag() > 0)
        {
          matchedElectronEta->Fill(scatEMomR.PseudoRapidity());
          electronEff->Fill(scatEMomR.PseudoRapidity());
          matchedElectronMomHist->Fill(scatEMomR.Mag());
          electronMomEff->Fill(scatEMomR.Mag());
          matchedElectronEtaMom->Fill(scatEMomR.PseudoRapidity(),scatEMomR.Mag());

          float deltaEta_electron = scatEMomT.PseudoRapidity() - scatEMomR.PseudoRapidity();
          float deltaPhi_electron = TVector2::Phi_mpi_pi(scatEMomT.Phi() - scatEMomR.Phi());
          float deltaR_electron = TMath::Sqrt(deltaEta_electron*deltaEta_electron + deltaPhi_electron*deltaPhi_electron);
          matchedElectronTrackDeltaR->Fill(deltaR_electron);

        }

        if (scatpMomR.Mag() > 0)
        {

          matchedProtonEta->Fill(scatpMomR.PseudoRapidity());
          protonEff->Fill(scatpMomR.PseudoRapidity());
          matchedProtonMomHist->Fill(scatpMomR.Mag());
          protonMomEff->Fill(scatpMomR.Mag());
          matchedProtonEtaMom->Fill(scatpMomR.PseudoRapidity(),scatpMomR.Mag());

          float deltaEta_proton = scatpMomT.PseudoRapidity() - scatpMomR.PseudoRapidity();
          float deltaPhi_proton = TVector2::Phi_mpi_pi(scatpMomT.Phi() - scatpMomR.Phi());
          float deltaR_proton = TMath::Sqrt(deltaEta_proton*deltaEta_proton + deltaPhi_proton*deltaPhi_proton);
          matchedProtonTrackDeltaR->Fill(deltaR_proton);
        }

        if (muPlusMomR.Mag() > 0)
        {
          
          matchedMuonEta->Fill(muPlusMomR.PseudoRapidity());
          muonEff->Fill(muPlusMomR.PseudoRapidity());
          matchedMuonMomHist->Fill(muPlusMomR.Mag());
          muonMomEff->Fill(muPlusMomR.Mag());
          matchedMuonEtaMom->Fill(muPlusMomR.PseudoRapidity(),muPlusMomR.Mag());

          float deltaEta_muPlus = muPlusMomT.PseudoRapidity() - muPlusMomR.PseudoRapidity();
          float deltaPhi_muPlus = TVector2::Phi_mpi_pi(muPlusMomT.Phi() - muPlusMomR.Phi());
          float deltaR_muPlus = TMath::Sqrt(deltaEta_muPlus*deltaEta_muPlus + deltaPhi_muPlus*deltaPhi_muPlus);
          matchedMuonTrackDeltaR->Fill(deltaR_muPlus);

        }

        if (muMinusMomR.Mag() > 0)
        {

          matchedMuonEta->Fill(muMinusMomR.PseudoRapidity());
          muonEff->Fill(muMinusMomR.PseudoRapidity());
          matchedMuonMomHist->Fill(muMinusMomR.Mag());
          muonMomEff->Fill(muMinusMomR.Mag());
          matchedMuonEtaMom->Fill(muMinusMomR.PseudoRapidity(),muMinusMomR.Mag()); 

          float deltaEta_muMinus = muMinusMomT.PseudoRapidity() - muMinusMomR.PseudoRapidity();
          float deltaPhi_muMinus = TVector2::Phi_mpi_pi(muMinusMomT.Phi() - muMinusMomR.Phi());
          float deltaR_muMinus = TMath::Sqrt(deltaEta_muMinus*deltaEta_muMinus + deltaPhi_muMinus*deltaPhi_muMinus);
          matchedMuonTrackDeltaR->Fill(deltaR_muMinus);
        }

        if (muPlusMomR.Mag() > 0 && muMinusMomR.Mag() > 0) // Check if the reconstructed muons have been found
        {
          JPsiMomR = muPlusMomR + muMinusMomR; 
          JPsi4MomR = muPlus4MomR + muMinus4MomR;
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

        }

        if (scatEMomR.Mag() > 0 && muPlusMomR.Mag() > 0 && muMinusMomR.Mag() > 0 && scatpMomR.Mag() > 0){ // if e', p', mu+, and mu- are in coincidence
  
          // Calculate kinematic variable

          double delta_h = (JPsi4MomR.E() + scatp4MomR.E()) - (JPsi4MomR.Pz() + scatp4MomR.Pz());
          double pt2_h = (pow((JPsi4MomR.Px()+scatp4MomR.Px()),2))+(pow((JPsi4MomR.Py()+scatp4MomR.Py()),2));
          double alpha_e = tan((scatE4MomR.Theta()/2));
          double alpha_h = delta_h/(sqrt(pt2_h));

          ROOT::Math::PxPyPzEVector vec_t_XPT_sum = (scatE4MomR + JPsi4MomR);
          TVector3 vec_t_XPT3 = TVector3(vec_t_XPT_sum.Px(), vec_t_XPT_sum.Py(), vec_t_XPT_sum.Pz());
          vec_t_XPT3.RotateY(0.025);
          ROOT::Math::PxPyPzEVector vec_t_XPT = ROOT::Math::PxPyPzEVector(vec_t_XPT3.X(), vec_t_XPT3.Y(), vec_t_XPT3.Z(), vec_t_XPT_sum.E());

          ROOT::Math::PxPyPzEVector vec_PMiss_Rec = (beamE4Mom + beamp4Mom) - (scatE4MomR + JPsi4MomR);
          ROOT::Math::PxPyPzEVector scatp4MomR_corr = ROOT::Math::PxPyPzEVector(vec_PMiss_Rec.P()*sin(scatp4MomT.Theta())*cos(scatp4MomT.Phi()), vec_PMiss_Rec.P()*sin(scatp4MomT.Theta())*sin(scatp4MomT.Phi()), vec_PMiss_Rec.P()*cos(scatp4MomT.Theta()), sqrt(pow(vec_PMiss_Rec.P(),2)+(pow(pMass,2))));

          y_e =(beamp4Mom.Dot(beamE4Mom - scatE4MomR))/(beamE4Mom.Dot(beamE4Mom));
          Q2_e = -(beamE4Mom - scatE4MomR).mag2();
          x_e = Q2_e/(4*beamE4Mom.E()*beamE4Mom.E()*y_e);

          y_JB = delta_h/(2*beamE4Mom.E());
          Q2_JB = pt2_h/(1-y_JB);
          x_JB = Q2_JB/(4*beamE4Mom.E()*beamE4Mom.E()*y_JB);

          y_DA = (alpha_h)/(alpha_e + alpha_h);
          Q2_DA = (4*beamE4Mom.E()*beamE4Mom.E())/(alpha_e*(alpha_e +alpha_h));
          x_DA = Q2_DA/(4*beamE4Mom.E()*beamE4Mom.E()*y_DA);
          
          y_sigma = delta_h/(delta_h + (scatE4MomR.E()*(1-cos(scatE4MomR.Theta()))));
          Q2_sigma = (pow(scatE4MomR.E(),2)*pow(sin(scatE4MomR.Theta()),2))/(1-y_sigma);
          x_sigma = Q2_sigma/(4*beamE4Mom.E()*beamE4Mom.E()*y_sigma);

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
          deltaQ2_e->Fill((Q2_e - Q2_truth)/Q2_truth);
          deltaQ2_JB->Fill((Q2_JB - Q2_truth)/Q2_truth);
          deltaQ2_DA->Fill((Q2_DA - Q2_truth)/Q2_truth);
          deltaQ2_sigma->Fill((Q2_sigma - Q2_truth)/Q2_truth);

          reconQ2_DA_vs_trueQ2->Fill(Q2_truth, Q2_DA);

          truet->Fill(t_truth);
          recont_eXBABE->Fill(t_eXBABE);
          recont_eXPT->Fill(t_eXPT);
          recont_eX->Fill(t_eX);
          recont_BABE->Fill(t_BABE);
          deltat_eXBABE->Fill((t_eXBABE - t_truth)/t_truth);
          deltat_eXPT->Fill((t_eXPT - t_truth)/t_truth);
          deltat_eX->Fill((t_eX - t_truth)/t_truth);
          deltat_BABE->Fill((t_BABE - t_truth)/t_truth);

          recont_eXBABE_vs_truet->Fill(t_truth, t_eXBABE);

          truex->Fill(x_truth);
          reconx_e->Fill(x_e);
          reconx_JB->Fill(x_JB);
          reconx_DA->Fill(x_DA);
          reconx_sigma->Fill(x_sigma);
          deltax_e->Fill((x_e - x_truth)/x_truth);
          deltax_JB->Fill((x_JB - x_truth)/x_truth);
          deltax_DA->Fill((x_DA - x_truth)/x_truth);
          deltax_sigma->Fill((x_sigma - x_truth)/x_truth);
        
        }

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

    }


    std::cout << "Event Processing Complete" << std::endl;
    
    
    electronEff->Divide(electronEta);
    protonEff->Divide(protonEta);
    muonEff->Divide(muonEta);
    JPsiEff->Divide(JPsiEta);

    electronMomEff->Divide(electronMomHist);
    protonMomEff->Divide(protonMomHist);
    muonMomEff->Divide(muonMomHist);
    JPsiMomEff->Divide(JPsiMomHist);

    truet_XbjkA->SetMarkerStyle(20);
    truet_XbjkB->SetMarkerStyle(20);
    truet_XbjkC->SetMarkerStyle(20);
    recont_XbjkA->SetMarkerStyle(20);
    recont_XbjkB->SetMarkerStyle(20);
    recont_XbjkC->SetMarkerStyle(20);

    truet_XbjkA->SetMarkerColor(kRed);
    truet_XbjkB->SetMarkerColor(kBlue);
    truet_XbjkC->SetMarkerColor(kGreen);
    recont_XbjkA->SetMarkerColor(kRed);
    recont_XbjkB->SetMarkerColor(kBlue);
    recont_XbjkC->SetMarkerColor(kGreen);

    truet_XbjkA->SetDrawOption("E1P");
    truet_XbjkB->SetDrawOption("E1P");
    truet_XbjkC->SetDrawOption("E1P");
    recont_XbjkA->SetDrawOption("E1P"); 
    recont_XbjkB->SetDrawOption("E1P");
    recont_XbjkC->SetDrawOption("E1P");

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
    matchedProtonEtaMom->Write();
    matchedProtonTrackDeltaR->Write();
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
    invMassMuMu->Write();
    matchedMuMuInvMass->Write();
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
    truex->Write();
    reconx_e->Write();
    reconx_JB->Write();
    reconx_DA->Write();
    reconx_sigma->Write();
    ofile->cd("..");
    ofile->mkdir("kinematicDifferences");
    ofile->cd("kinematicsDifferences");
    deltaQ2_e->Write();
    deltaQ2_JB->Write();
    deltaQ2_DA->Write();
    deltaQ2_sigma->Write();
    deltat_eXBABE->Write();
    deltat_eXPT->Write();
    deltat_eX->Write();
    deltat_BABE->Write();
    deltax_e->Write();
    deltax_JB->Write();
    deltax_DA->Write();
    deltax_sigma->Write();
    ofile->cd("..");

    ofile->Close(); // Close output file

    std::cout << "Histograms written to file: " << outfilename << std::endl;
    std::cout << "Track analysis completed." << std::endl;
  }