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
#include <TPaveText.h>
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
#include <TLegend.h>
#include <TLine.h>

#include "ePICStyle.C"

void DVMP_histograms()
{
    //gROOT->SetBatch(kTRUE);
    gROOT->ProcessLine("SetePICStyle()");
    gStyle->SetOptStat(0);

    gStyle->SetCanvasPreferGL(true);


    std::vector<int> sixColourScheme = {
    TColor::GetColor("#7021dd"),     // violet
    TColor::GetColor("#964a8b"),     // grape
    TColor::GetColor("#e42536"),     // red
    TColor::GetColor("#f89c20"),     // yellow
    TColor::GetColor("#5790fc"),     // blue
    TColor::GetColor("#9c9ca1"),     // grey
    };

    TString infileN = "outputs/DVMP_JPsi_AnalysisOutput_10ifb_10x130ep.root";
    TFile *inFile = new TFile(infileN);

    std::string outfilename = "outputs/DVMP_JPsi_hists_10ifb_10x130ep.root";

    // Set output file for the histograms
    TFile *ofile = TFile::Open(outfilename.c_str(),"RECREATE");

    // Load Invariant Mass, Missing Mass, and Angle Histograms

    TH1D *muonScatterAngle = (TH1D*) inFile->Get("invariantMass/muonScatterAngle");

    TH1D *invMassMuMu = (TH1D*) inFile->Get("invariantMass/invMassMuMu");
    TH1D *matchedMuMuInvMass = (TH1D*) inFile->Get("invariantMass/matchedMuMuInvMass");

    TH1D *matchedMissingMassEpToJPsiX = (TH1D*) inFile->Get("invariantMass/matchedMissingMassEpToJPsiX");
    TH1D *matchedMissingMass2EpToJPsiX = (TH1D*) inFile->Get("invariantMass/matchedMissingMass2EpToJPsiX");

    // Load Kinematics Histograms
    TH1D* trueQ2 = (TH1D*) inFile->Get("kinematics/trueQ2");
    TH1D* reconQ2_DA = (TH1D*) inFile->Get("kinematics/reconQ2_DA");
    TH1D* reconQ2_JB = (TH1D*) inFile->Get("kinematics/reconQ2_JB");
    TH1D* reconQ2_e = (TH1D*) inFile->Get("kinematics/reconQ2_e");
    TH1D* reconQ2_sigma = (TH1D*) inFile->Get("kinematics/reconQ2_sigma");

    TH1D* deltaQ2_DA = (TH1D*) inFile->Get("kinematicDifferences/deltaQ2_DA");
    TH1D* deltaQ2_JB = (TH1D*) inFile->Get("kinematicDifferences/deltaQ2_JB");
    TH1D* deltaQ2_e = (TH1D*) inFile->Get("kinematicDifferences/deltaQ2_e");
    TH1D* deltaQ2_sigma = (TH1D*) inFile->Get("kinematicDifferences/deltaQ2_sigma");

    TH2D* reconQ2_DA_vs_trueQ2 = (TH2D*) inFile->Get("kinematicDifferences/reconQ2_DA_vs_trueQ2");

    TH1D* truet = (TH1D*) inFile->Get("kinematics/truet");
    TH1D* recont_eXBABE = (TH1D*) inFile->Get("kinematics/recont_eXBABE");
    TH1D* recont_eXPT = (TH1D*) inFile->Get("kinematics/recont_eXPT");
    TH1D* recont_eX = (TH1D*) inFile->Get("kinematics/recont_eX");
    TH1D* recont_BABE = (TH1D*) inFile->Get("kinematics/recont_BABE");

    TH2D* recont_eXBABE_vs_truet = (TH2D*) inFile->Get("kinematics/recont_eXBABE_vs_truet");

    TH1D* deltat_eXBABE = (TH1D*) inFile->Get("kinematicDifferences/deltat_eXBABE");
    TH1D* deltat_eXPT = (TH1D*) inFile->Get("kinematicDifferences/deltat_eXPT");
    TH1D* deltat_eX = (TH1D*) inFile->Get("kinematicDifferences/deltat_eX");
    TH1D* deltat_BABE = (TH1D*) inFile->Get("kinematicDifferences/deltat_BABE");

    TH2D* deltat_eXBABE_vs_truet = (TH2D*) inFile->Get("kinematicDifferences/deltat_eXBABE_vs_truet");

    TH1D* truey = (TH1D*) inFile->Get("kinematics/truey");
    TH1D* recony_DA = (TH1D*) inFile->Get("kinematics/recony_DA");
    TH1D* recony_JB = (TH1D*) inFile->Get("kinematics/recony_JB");
    TH1D* recony_e = (TH1D*) inFile->Get("kinematics/recony_e");
    TH1D* recony_sigma = (TH1D*) inFile->Get("kinematics/recony_sigma");

    TH2D* recony_DA_vs_truey = (TH2D*) inFile->Get("kinematicDifferences/recony_DA_vs_truey");

    TH1D* deltay_DA = (TH1D*) inFile->Get("kinematicDifferences/deltay_DA");
    TH1D* deltay_JB = (TH1D*) inFile->Get("kinematicDifferences/deltay_JB");
    TH1D* deltay_e = (TH1D*) inFile->Get("kinematicDifferences/deltay_e");
    TH1D* deltay_sigma = (TH1D*) inFile->Get("kinematicDifferences/deltay_sigma");

    TH1D* truex = (TH1D*) inFile->Get("kinematics/truex");
    TH1D* reconx_DA = (TH1D*) inFile->Get("kinematics/reconx_DA");
    TH1D* reconx_JB = (TH1D*) inFile->Get("kinematics/reconx_JB");
    TH1D* reconx_e = (TH1D*) inFile->Get("kinematics/reconx_e");
    TH1D* reconx_sigma = (TH1D*) inFile->Get("kinematics/reconx_sigma");

    TH1D* deltax_DA = (TH1D*) inFile->Get("kinematicDifferences/deltax_DA");
    TH1D* deltax_JB = (TH1D*) inFile->Get("kinematicDifferences/deltax_JB");
    TH1D* deltax_e = (TH1D*) inFile->Get("kinematicDifferences/deltax_e");
    TH1D* deltax_sigma = (TH1D*) inFile->Get("kinematicDifferences/deltax_sigma");

    TH1D* truet_XbjkA = (TH1D*) inFile->Get("tDistributionsByXbjk/truet_XbjkA");
    TH1D* recont_XbjkA = (TH1D*) inFile->Get("tDistributionsByXbjk/recont_XbjkA");
    TH1D* truet_XbjkB = (TH1D*) inFile->Get("tDistributionsByXbjk/truet_XbjkB");
    TH1D* recont_XbjkB = (TH1D*) inFile->Get("tDistributionsByXbjk/recont_XbjkB");
    TH1D* truet_XbjkC = (TH1D*) inFile->Get("tDistributionsByXbjk/truet_XbjkC");
    TH1D* recont_XbjkC = (TH1D*) inFile->Get("tDistributionsByXbjk/recont_XbjkC");
    
    
    //
    // Muon angle
    //

    TCanvas* C_MuonAngle = new TCanvas("C_MuonAngle", "C_MuonAngle",1100,800);
    TLegend* legend_MuonAngle = new TLegend(0.55, 0.7, 0.8, 0.9);
    muonScatterAngle->Draw();
    legend_MuonAngle->AddEntry(muonScatterAngle, "Scattered angle", "l");

    TPaveText* text_MuonAngle = new TPaveText(0.2, 0.7, 0.4, 0.9, "NDC");
    text_MuonAngle->AddText("ePIC Performance");
    text_MuonAngle->AddText("e+p, 10ex130p");
    text_MuonAngle->AddText("lAger 3.6.1");
    text_MuonAngle->AddText("DVMP J/#psi -> #mu#mu");
    text_MuonAngle->SetTextAlign(12);
    text_MuonAngle->SetFillStyle(0);
    text_MuonAngle->SetBorderSize(0);

    TLine *MuonAngleL_line = new TLine(2.5, 0, 2.5, muonScatterAngle->GetBinContent(muonScatterAngle->GetMaximumBin()));
    MuonAngleL_line->SetLineColor(kRed);
    MuonAngleL_line->SetLineStyle(2);   // dashed (optional)
    MuonAngleL_line->SetLineWidth(2);
    MuonAngleL_line->Draw();

    TLine *MuonAngleU_line = new TLine(4, 0, 4, muonScatterAngle->GetBinContent(muonScatterAngle->GetMaximumBin()));
    MuonAngleU_line->SetLineColor(kRed);
    MuonAngleU_line->SetLineStyle(2);   // dashed (optional)
    MuonAngleU_line->SetLineWidth(2);
    MuonAngleU_line->Draw();


    C_MuonAngle->Update();
    legend_MuonAngle->Draw();
    text_MuonAngle->Draw();
    
    

    //
    // Missing mass
    //

    TCanvas* C_MissingMass = new TCanvas("C_MissingMass", "C_MissingMass",1100,800);
    TLegend* legend_MissingMass = new TLegend(0.6, 0.7, 0.8, 0.9);
    matchedMissingMassEpToJPsiX->Draw();
    legend_MissingMass->AddEntry(matchedMissingMassEpToJPsiX, "Missing Mass", "l");

    TPaveText* text_MissingMass = new TPaveText(0.2, 0.7, 0.4, 0.9, "NDC");
    text_MissingMass->AddText("ePIC Performance");
    text_MissingMass->AddText("e+p, 10ex130p");
    text_MissingMass->AddText("lAger 3.6.1");
    text_MissingMass->AddText("DVMP J/#psi -> #mu#mu");
    text_MissingMass->SetTextAlign(12);
    text_MissingMass->SetFillStyle(0);
    text_MissingMass->SetBorderSize(0);

    TLine *MM_line = new TLine(2, 0, 2, matchedMissingMassEpToJPsiX->GetBinContent(matchedMissingMassEpToJPsiX->GetMaximumBin()));
    MM_line->SetLineColor(kRed);
    MM_line->SetLineStyle(2);   // dashed (optional)
    MM_line->SetLineWidth(2);
    MM_line->Draw();

    C_MissingMass->SetLogy();
    C_MissingMass->Update();
    legend_MissingMass->Draw();
    text_MissingMass->Draw();

    //
    // Missing mass squared
    //

    TCanvas* C_MissingMass2 = new TCanvas("C_MissingMass2", "C_MissingMass2",1100,800);
    TLegend* legend_MissingMass2 = new TLegend(0.5, 0.7, 0.7, 0.9);
    matchedMissingMass2EpToJPsiX->Draw();
    legend_MissingMass2->AddEntry(matchedMissingMass2EpToJPsiX, "Missing Mass Squared", "l");

    TPaveText* text_MissingMass2 = new TPaveText(0.2, 0.7, 0.4, 0.9, "NDC");
    text_MissingMass2->AddText("ePIC Performance");
    text_MissingMass2->AddText("e+p, 10ex130p");
    text_MissingMass2->AddText("lAger 3.6.1");
    text_MissingMass2->AddText("DVMP J/#psi -> #mu#mu");
    text_MissingMass2->SetTextAlign(12);
    text_MissingMass2->SetFillStyle(0);
    text_MissingMass2->SetBorderSize(0);

    TLine *MM2_line = new TLine(-1, 0, -1, matchedMissingMass2EpToJPsiX->GetBinContent(matchedMissingMass2EpToJPsiX->GetMaximumBin()));
    MM2_line->SetLineColor(kRed);
    MM2_line->SetLineStyle(2);   // dashed (optional)
    MM2_line->SetLineWidth(2);
    MM2_line->Draw();

    C_MissingMass2->SetLogy();
    C_MissingMass2->Update();
    legend_MissingMass2->Draw();
    text_MissingMass2->Draw();

    //
    // True Invariant mass
    //

    TCanvas* C_InvMassTrue = new TCanvas("C_InvMassTrue", "C_InvMassTrue",1100,800);
    TLegend* legend_InvMassTrue = new TLegend(0.2, 0.3, 0.4, 0.6);
    invMassMuMu->Draw();
    legend_InvMassTrue->AddEntry(invMassMuMu, "Invariant Mass of #mu#mu", "l");

    TPaveText* text_InvMassTrue = new TPaveText(0.2, 0.7, 0.4, 0.9, "NDC");
    text_InvMassTrue->AddText("ePIC Performance");
    text_InvMassTrue->AddText("e+p, 10ex130p");
    text_InvMassTrue->AddText("lAger 3.6.1");
    text_InvMassTrue->AddText("DVMP J/#psi -> #mu#mu");
    text_InvMassTrue->SetTextAlign(12);
    text_InvMassTrue->SetFillStyle(0);
    text_InvMassTrue->SetBorderSize(0);

    C_InvMassTrue->Update();
    legend_InvMassTrue->Draw();
    text_InvMassTrue->Draw();
    
    
    //
    // Rec Invariant mass
    //

    TCanvas* C_InvMassRec = new TCanvas("C_InvMassRec", "C_InvMassRec",1100,800);
    TLegend* legend_InvMassRec = new TLegend(0.2, 0.3, 0.4, 0.6);
    matchedMuMuInvMass->Draw();
    legend_InvMassRec->AddEntry(matchedMuMuInvMass, "Invariant Mass of #mu#mu", "l");

    TPaveText* text_InvMassRec = new TPaveText(0.2, 0.7, 0.4, 0.9, "NDC");
    text_InvMassRec->AddText("ePIC Performance");
    text_InvMassRec->AddText("e+p, 10ex130p");
    text_InvMassRec->AddText("lAger 3.6.1");
    text_InvMassRec->AddText("DVMP J/#psi -> #mu#mu");
    text_InvMassRec->SetTextAlign(12);
    text_InvMassRec->SetFillStyle(0);
    text_InvMassRec->SetBorderSize(0);

    TLine *InvMassL_line = new TLine(2.5, 0, 2.5, matchedMuMuInvMass->GetBinContent(matchedMuMuInvMass->GetMaximumBin()));
    InvMassL_line->SetLineColor(kRed);
    InvMassL_line->SetLineStyle(2);   // dashed (optional)
    InvMassL_line->SetLineWidth(2);
    InvMassL_line->Draw();

    TLine *InvMassU_line = new TLine(4, 0, 4, matchedMuMuInvMass->GetBinContent(matchedMuMuInvMass->GetMaximumBin()));
    InvMassU_line->SetLineColor(kRed);
    InvMassU_line->SetLineStyle(2);   // dashed (optional)
    InvMassU_line->SetLineWidth(2);
    InvMassU_line->Draw();


    C_InvMassRec->Update();
    legend_InvMassRec->Draw();
    text_InvMassRec->Draw();
    
    
    //
    // Q2 diff
    //

    TCanvas* C_Q2diff = new TCanvas("C_Q2diff", "C_Q2diff",1100,800);
    TLegend* legend_Q2diff = new TLegend(0.5, 0.7, 0.7, 0.9);

    deltaQ2_DA->SetFillColorAlpha(sixColourScheme[0], 0.75);
    deltaQ2_DA->Draw("HIST F");
    legend_Q2diff->AddEntry(deltaQ2_DA, "DA Method", "f");
    deltaQ2_JB->SetFillColorAlpha(sixColourScheme[2], 0.75);
    deltaQ2_JB->Draw("HIST F SAME");
    legend_Q2diff->AddEntry(deltaQ2_JB, "JB Method", "f");
    deltaQ2_e->SetFillColorAlpha(sixColourScheme[3], 0.75);
    deltaQ2_e->Draw("HIST F SAME");
    legend_Q2diff->AddEntry(deltaQ2_e, "Electron Method", "f");
    deltaQ2_sigma->SetFillColorAlpha(sixColourScheme[4], 0.75);
    deltaQ2_sigma->Draw("HIST F SAME");
    legend_Q2diff->AddEntry(deltaQ2_sigma, "Sigma Method", "f");

    TPaveText* text_Q2diff = new TPaveText(0.2, 0.7, 0.4, 0.9, "NDC");
    text_Q2diff->AddText("ePIC Performance");
    text_Q2diff->AddText("e+p, 10ex130p");
    text_Q2diff->AddText("lAger 3.6.1");
    text_Q2diff->AddText("DVMP J/#psi -> #mu#mu");
    text_Q2diff->SetTextAlign(12);
    text_Q2diff->SetFillStyle(0);
    text_Q2diff->SetBorderSize(0);

    C_Q2diff->Update();
    legend_Q2diff->Draw();
    text_Q2diff->Draw();

    
    //
    // t diff
    //


    TCanvas* C_tdiff = new TCanvas("C_tdiff", "C_tdiff",1100,800);
    TLegend* legend_tdiff = new TLegend(0.5, 0.7, 0.7, 0.9);

    deltat_eXBABE->SetFillColorAlpha(sixColourScheme[0], 0.75);
    deltat_eXBABE->Draw("HIST F");
    legend_tdiff->AddEntry(deltat_eXBABE, "eXBABE Method", "f");
    deltat_eXPT->SetFillColorAlpha(sixColourScheme[2], 0.75);
    deltat_eXPT->Draw("HIST F SAME");
    legend_tdiff->AddEntry(deltat_eXPT, "eXPT Method", "f");
    deltat_eX->SetFillColorAlpha(sixColourScheme[3], 0.75);
    deltat_eX->Draw("HIST F SAME");
    legend_tdiff->AddEntry(deltat_eX, "eX Method", "f");
    deltat_BABE->SetFillColorAlpha(sixColourScheme[4], 0.75);
    deltat_BABE->Draw("HIST F SAME");
    legend_tdiff->AddEntry(deltat_BABE, "BABE Method", "f");

    TPaveText* text_tdiff = new TPaveText(0.2, 0.7, 0.4, 0.9, "NDC");
    text_tdiff->AddText("ePIC Performance");
    text_tdiff->AddText("e+p, 10ex130p");
    text_tdiff->AddText("lAger 3.6.1");
    text_tdiff->AddText("DVMP J/#psi -> #mu#mu");
    text_tdiff->SetTextAlign(12);
    text_tdiff->SetFillStyle(0);
    text_tdiff->SetBorderSize(0);

    C_tdiff->Update();
    legend_tdiff->Draw();
    text_tdiff->Draw();

    
    //
    // t_xbjk
    //

    TCanvas* C_t_xbjk = new TCanvas("C_t_xbjk", "C_t_xbjk",1100,800);
    TLegend* legend_t_xbjk = new TLegend(0.5, 0.7, 0.7, 0.9);

    recont_XbjkA->SetMarkerColor(sixColourScheme[0]);
    recont_XbjkA->Draw("ELP");
    legend_t_xbjk->AddEntry(recont_XbjkA, "0.0016 < x_{bjk} < 0.0025", "p");
    recont_XbjkB->SetMarkerColor(sixColourScheme[2]);
    recont_XbjkB->Draw("ELP SAME");
    legend_t_xbjk->AddEntry(recont_XbjkB, "0.016 < x_{bjk} < 0.025", "p");
    recont_XbjkC->SetMarkerColor(sixColourScheme[4]);
    recont_XbjkC->Draw("ELP SAME");
    legend_t_xbjk->AddEntry(recont_XbjkC, "0.16 < x_{bjk} < 0.25", "p");

    TPaveText* text_t_xbjk = new TPaveText(0.575, 0.495, 0.775, 0.695, "NDC");
    text_t_xbjk->AddText("ePIC Performance");
    text_t_xbjk->AddText("e+p, 10ex130p");
    text_t_xbjk->AddText("lAger 3.6.1");
    text_t_xbjk->AddText("DVMP J/#psi -> #mu#mu");
    text_t_xbjk->SetTextAlign(12);
    text_t_xbjk->SetFillStyle(0);
    text_t_xbjk->SetBorderSize(0);

    C_t_xbjk->SetLogy();
    C_t_xbjk->Update();
    legend_t_xbjk->Draw();
    text_t_xbjk->Draw();


}