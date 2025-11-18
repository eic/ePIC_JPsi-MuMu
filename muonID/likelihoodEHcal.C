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

#include "../ePICStyle.C"

void likelihoodEHcal()
{
    //gROOT->SetBatch(kTRUE);
    //gROOT->ProcessLine("SetePICStyle()");
    gStyle->SetOptStat(0);
    gStyle->SetCanvasPreferGL(true); // enables transparency

    std::vector<int> sixColourScheme = {
    TColor::GetColor("#7021dd"),     // violet
    TColor::GetColor("#964a8b"),     // grape
    TColor::GetColor("#e42536"),     // red
    TColor::GetColor("#f89c20"),     // yellow
    TColor::GetColor("#5790fc"),     // blue
    TColor::GetColor("#9c9ca1"),     // grey
    };

    TString infileN = "outputs/mu-pi_epPlots.root";
    TFile *inFile = new TFile(infileN);

    std::string outfilename = "outputs/mu-pi_likelihood.root";

    // Set output file for the histograms
    TFile *ofile = TFile::Open(outfilename.c_str(),"RECREATE");

    // Ecal ep Plots by particle
    TH2D *pionEpTrueEcal = (TH2D*) inFile->Get("trueEcalEp/pionEpTrueEcal");
    TH2D *muonEpTrueEcal = (TH2D*) inFile->Get("trueEcalEp/muonEpTrueEcal");

    TH2D *pionEpRecEcal = (TH2D*) inFile->Get("recoEcalEp/pionEpRecEcal");
    TH2D *muonEpRecEcal = (TH2D*) inFile->Get("recoEcalEp/muonEpRecEcal");

    // Hcal ep Plots by particle
    TH2D *pionEpTrueHcal = (TH2D*) inFile->Get("trueHcalEp/pionEpTrueHcal");
    TH2D *muonEpTrueHcal = (TH2D*) inFile->Get("trueHcalEp/muonEpTrueHcal");

    TH2D *pionEpRecHcal = (TH2D*) inFile->Get("recoHcalEp/protonEpRecHcal");
    TH2D *muonEpRecHcal = (TH2D*) inFile->Get("recoHcalEp/protonEpRecHcal");


    TH1D *EpTrueEcal[20][2];
    TH1D *EpRecEcal[20][2];
    TH1D *EpTrueHcal[20][2];
    TH1D *EpRecHcal[20][2];

    std::string particleNameList[2] = {"Muon", "Pion"};

    for (int particleCount = 0; particleCount < 2; particleCount++)
    {

        for (int momCount = 0; momCount < 20; momCount++)
        {
            std::string histN = particleNameList[particleCount] + std::string("EpTrueEcalMomBin") + std::to_string(momCount);
            TString histNT = histN;
            std::string histT = particleNameList[particleCount] + std::string("True Ecal E/p for Momentum Bin: ") + std::to_string(std::round((momCount+0.5)*10)/10) + " GeV";
            TString histTT = histT;

            EpTrueEcal[momCount][particleCount] = new TH1D(histNT, histTT, 100.,0.,2.0);
        }

        for (int momCount = 0; momCount < 20; momCount++)
        {
            std::string histN = particleNameList[particleCount] + std::string("EpRecEcalMomBin") + std::to_string(momCount);
            TString histNT = histN;
            std::string histT = particleNameList[particleCount] + std::string("Reconstructed Ecal E/p for Momentum Bin: ") + std::to_string(std::round((momCount+0.5)*10)/10) + " GeV";
            TString histTT = histT;

            EpRecEcal[momCount][particleCount] = new TH1D(histNT, histTT, 100.,0.,2.0);
        }

        for (int momCount = 0; momCount < 20; momCount++)
        {
            std::string histN = particleNameList[particleCount] + std::string("EpTrueHcalMomBin") + std::to_string(momCount);
            TString histNT = histN;
            std::string histT = particleNameList[particleCount] + std::string("True Hcal E/p for Momentum Bin: ") + std::to_string(std::round((momCount+0.5)*10)/10) + " GeV";
            TString histTT = histT;

            EpTrueHcal[momCount][particleCount] = new TH1D(histNT, histTT, 100.,0.,2.0);
        }

        for (int momCount = 0; momCount < 20; momCount++)
        {
            std::string histN = particleNameList[particleCount] + std::string("EpRecHcalMomBin") + std::to_string(momCount);
            TString histNT = histN;
            std::string histT = particleNameList[particleCount] + std::string("Reconstructed Hcal E/p for Momentum Bin: ") + std::to_string(std::round((momCount+0.5)*10)/10) + " GeV";
            TString histTT = histT;

            EpRecHcal[momCount][particleCount] = new TH1D(histNT, histTT, 100.,0.,2.0);
        }

    }

    for (int binCount = 1; binCount < muonEpTrueEcal->GetNbinsX(); binCount++)
    {

        int momBin = int(muonEpTrueEcal->GetXaxis()->GetBinCenter(binCount));

        if (momBin >= 20) break;

        for (int yBinCount = 1; yBinCount < muonEpTrueEcal->GetNbinsY(); yBinCount++)
        {
            EpTrueEcal[momBin][0]->Fill(muonEpTrueEcal->GetYaxis()->GetBinCenter(yBinCount), muonEpTrueEcal->GetBinContent(binCount, yBinCount));
            EpTrueEcal[momBin][1]->Fill(pionEpTrueEcal->GetYaxis()->GetBinCenter(yBinCount), pionEpTrueEcal->GetBinContent(binCount, yBinCount));

            EpRecEcal[momBin][0]->Fill(muonEpRecEcal->GetYaxis()->GetBinCenter(yBinCount), muonEpRecEcal->GetBinContent(binCount, yBinCount));
            EpRecEcal[momBin][1]->Fill(pionEpRecEcal->GetYaxis()->GetBinCenter(yBinCount), pionEpRecEcal->GetBinContent(binCount, yBinCount));

            EpTrueHcal[momBin][0]->Fill(muonEpTrueHcal->GetYaxis()->GetBinCenter(yBinCount), muonEpTrueHcal->GetBinContent(binCount, yBinCount));
            EpTrueHcal[momBin][1]->Fill(pionEpTrueHcal->GetYaxis()->GetBinCenter(yBinCount), pionEpTrueHcal->GetBinContent(binCount, yBinCount));

            EpRecHcal[momBin][0]->Fill(muonEpRecHcal->GetYaxis()->GetBinCenter(yBinCount), muonEpRecHcal->GetBinContent(binCount, yBinCount));
            EpRecHcal[momBin][1]->Fill(pionEpRecHcal->GetYaxis()->GetBinCenter(yBinCount), pionEpRecHcal->GetBinContent(binCount, yBinCount));

        }
    }

    TCanvas *canvas_EpTrueEcal = new TCanvas("canvas_EpTrueEcal","canvas_EpTrueEcal",1100,800);
    canvas_EpTrueEcal->Divide(4,5);
    for (int momCount = 0; momCount < 20; momCount++)
    {
        EpTrueEcal[momCount][0]->SetLineColor(sixColourScheme[4]);
        EpTrueEcal[momCount][0]->SetFillStyle(1001);
        EpTrueEcal[momCount][0]->SetFillColorAlpha(sixColourScheme[4], 0.25);

        EpTrueEcal[momCount][1]->SetLineColor(sixColourScheme[2]);
        EpTrueEcal[momCount][1]->SetFillStyle(1001);
        EpTrueEcal[momCount][1]->SetFillColorAlpha(sixColourScheme[2], 0.25);

        canvas_EpTrueEcal->cd(momCount+1);
        EpTrueEcal[momCount][0]->Draw("HIST");
        EpTrueEcal[momCount][1]->Draw("HIST SAME");
    }

    canvas_EpTrueEcal->Update();

    TCanvas *canvas_EpRecEcal = new TCanvas("canvas_EpRecEcal","canvas_EpRecEcal",1100,800);
    canvas_EpRecEcal->Divide(4,5);
    for (int momCount = 0; momCount < 20; momCount++)
    {
        EpRecEcal[momCount][0]->SetLineColor(sixColourScheme[4]);
        EpRecEcal[momCount][0]->SetFillStyle(1001);
        EpRecEcal[momCount][0]->SetFillColorAlpha(sixColourScheme[4], 0.25);

        EpRecEcal[momCount][1]->SetLineColor(sixColourScheme[2]);
        EpRecEcal[momCount][1]->SetFillStyle(1001);
        EpRecEcal[momCount][1]->SetFillColorAlpha(sixColourScheme[2], 0.25);

        canvas_EpRecEcal->cd(momCount+1);
        EpRecEcal[momCount][0]->Draw("HIST");
        EpRecEcal[momCount][1]->Draw("HIST SAME");
    }

    canvas_EpRecEcal->Update();

    TCanvas *canvas_EpTrueHcal = new TCanvas("canvas_EpTrueHcal","canvas_EpTrueHcal",1100,800);
    canvas_EpTrueHcal->Divide(4,5);
    for (int momCount = 0; momCount < 20; momCount++)
    {
        EpTrueHcal[momCount][0]->SetLineColor(sixColourScheme[4]);
        EpTrueHcal[momCount][0]->SetFillStyle(1001);
        EpTrueHcal[momCount][0]->SetFillColorAlpha(sixColourScheme[4], 0.25);

        EpTrueHcal[momCount][1]->SetLineColor(sixColourScheme[2]);
        EpTrueHcal[momCount][1]->SetFillStyle(1001);
        EpTrueHcal[momCount][1]->SetFillColorAlpha(sixColourScheme[2], 0.25);

        canvas_EpTrueHcal->cd(momCount+1);
        EpTrueHcal[momCount][0]->Draw("HIST");
        EpTrueHcal[momCount][1]->Draw("HIST SAME");
    }

    canvas_EpTrueHcal->Update();

    TCanvas *canvas_EpRecHcal = new TCanvas("canvas_EpRecHcal","canvas_EpRecHcal",1100,800);
    canvas_EpRecHcal->Divide(4,5);
    for (int momCount = 0; momCount < 20; momCount++)
    {
        EpRecHcal[momCount][0]->SetLineColor(sixColourScheme[4]);
        EpRecHcal[momCount][0]->SetFillStyle(1001);
        EpRecHcal[momCount][0]->SetFillColorAlpha(sixColourScheme[4], 0.25);

        EpRecHcal[momCount][1]->SetLineColor(sixColourScheme[2]);
        EpRecHcal[momCount][1]->SetFillStyle(1001);
        EpRecHcal[momCount][1]->SetFillColorAlpha(sixColourScheme[2], 0.25);

        canvas_EpRecHcal->cd(momCount+1);
        EpRecHcal[momCount][0]->Draw("HIST");
        EpRecHcal[momCount][1]->Draw("HIST SAME");
    }

    canvas_EpRecHcal->Update();


    ofile->cd();
    ofile->mkdir("EpTrueEcal");
    ofile->cd("EpTrueEcal");
    for (int momCount = 0; momCount < 20; momCount++)
    {
        EpTrueEcal[momCount][0]->Write();
        EpTrueEcal[momCount][1]->Write();
    }
    ofile->cd("..");
    ofile->mkdir("EpRecEcal");
    ofile->cd("EpRecEcal");
    for (int momCount = 0; momCount < 20; momCount++)
    {
        EpRecEcal[momCount][0]->Write();
        EpRecEcal[momCount][1]->Write();
    }
    ofile->cd("..");
    ofile->mkdir("EpTrueHcal");
    ofile->cd("EpTrueHcal");
    for (int momCount = 0; momCount < 20; momCount++)
    {
        EpTrueHcal[momCount][0]->Write();
        EpTrueHcal[momCount][1]->Write();
    }
    ofile->cd("..");
    ofile->mkdir("EpRecHcal");
    ofile->cd("EpRecHcal");
    for (int momCount = 0; momCount < 20; momCount++)
    {
        EpRecHcal[momCount][0]->Write();
        EpRecHcal[momCount][1]->Write();
    }

    ofile->cd("..");
    ofile->Close();
}