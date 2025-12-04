#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TVector3.h>
#include <TMath.h>
#include <TF1.h>

#include "DVMP_JPsi_Analysis.h"


bool IsMuon(TVector3 trackMom, int simuID,
            TTreeReaderArray<float>& EcalBarrelEng, TTreeReaderArray<float>& EcalEndcapPEng, TTreeReaderArray<float>& EcalEndcapNEng, 
            TTreeReaderArray<float>& HcalBarrelEng, TTreeReaderArray<float>& HcalEndcapPEng, TTreeReaderArray<float>& LFHcalEng, TTreeReaderArray<float>& HcalEndcapNEng, 
            TTreeReaderArray<unsigned int>& simuAssocEcalBarrel, TTreeReaderArray<unsigned int>& simuAssocEcalEndcapP, TTreeReaderArray<unsigned int>& simuAssocEcalEndcapN,
            TTreeReaderArray<unsigned int>& simuAssocHcalBarrel, TTreeReaderArray<unsigned int>& simuAssocHcalEndcapP, TTreeReaderArray<unsigned int>& simuAssocLFHcal, TTreeReaderArray<unsigned int>& simuAssocHcalEndcapN)
{

    double ECalEnergy = 0.0;
    double HCalEnergy = 0.0;
    int ECalHits = 0;
    int HCalHits = 0;

    // E/p cut functions

    TF1 *eCalCutFunction = new TF1("eCalCutFunction","[0]/(x + [1])",0.,20.);
    eCalCutFunction->SetParameter(0,0.29);
    eCalCutFunction->SetParameter(1,0.23);

    TF1 *eCalCutWidthFunction = new TF1("eCalCutWidthFunction","3*([0]/x)",0.,20.);
    eCalCutWidthFunction->SetParameter(0,0.09);

    TF1 *hCalCutFunction = new TF1("hCalCutFunction","[0]/(x + [1])",0.,20.);
    hCalCutFunction->SetParameter(0,1.03);
    hCalCutFunction->SetParameter(1,0.44);

    TF1 *hCalCutWidthFunction = new TF1("hCalCutWidthFunction","3*([0]/x)",0.,20.);
    hCalCutWidthFunction->SetParameter(0,0.21);

    if (EcalBarrelEng.GetSize() == simuAssocEcalBarrel.GetSize())
    {
        for (int jB = 0; jB < simuAssocEcalBarrel.GetSize(); jB++) // Look for associations in the Ecal Barrel
        {
            if (simuAssocEcalBarrel[jB] == simuID)
            {
                ECalEnergy += EcalBarrelEng[jB];
                ECalHits += 1.;
            }
        }
    }
    if (EcalEndcapPEng.GetSize() == simuAssocEcalEndcapP.GetSize()) 
    {
        for (int jP = 0; jP < simuAssocEcalEndcapP.GetSize(); jP++) // Look for associations in the Ecal Endcap P
        {
            if (simuAssocEcalEndcapP[jP] == simuID)
            {
                ECalEnergy += EcalEndcapPEng[jP];
                ECalHits += 1.;
            }
        }
    }
    if (EcalEndcapNEng.GetSize() == simuAssocEcalEndcapN.GetSize())
    {
        for (int jN = 0; jN < simuAssocEcalEndcapN.GetSize(); jN++) // Look for associations in the Ecal Endcap N
        {
            if (simuAssocEcalEndcapN[jN] == simuID)
            {
                ECalEnergy += EcalEndcapNEng[jN];
                ECalHits += 1.;
            }
        }
    }

    if (HcalBarrelEng.GetSize() == simuAssocHcalBarrel.GetSize())
    {
        for (int jB = 0; jB < simuAssocHcalBarrel.GetSize(); jB++) // Look for associations in the Ecal Barrel
        {
            if (simuAssocHcalBarrel[jB] == simuID)
            {
                HCalEnergy += HcalBarrelEng[jB];
                HCalHits += 1.;
            }
        }
    }
            
    if (HcalEndcapPEng.GetSize() == simuAssocHcalEndcapP.GetSize())
    {
        for (int jP = 0; jP < simuAssocHcalEndcapP.GetSize(); jP++) // Look for associations in the Ecal Endcap P
        {
            if (simuAssocHcalEndcapP[jP] == simuID)
            {
                HCalEnergy += HcalEndcapPEng[jP];
                HCalHits += 1.;
            }
        }
    }

    if (LFHcalEng.GetSize() == simuAssocLFHcal.GetSize())
    {
        for (int jL = 0; jL < simuAssocLFHcal.GetSize(); jL++) // Look for associations in the LFHcal
        {
            if (simuAssocLFHcal[jL] == simuID)
            {
                HCalEnergy += LFHcalEng[jL];
                HCalHits += 1.;
            }
        }
    }

            
    if (HcalEndcapNEng.GetSize() == simuAssocHcalEndcapN.GetSize())
    {
        for (int jN = 0; jN < simuAssocHcalEndcapN.GetSize(); jN++) // Look for associations in the Ecal Endcap N
        {
            if (simuAssocHcalEndcapN[jN] == simuID)
            {
                HCalEnergy += HcalEndcapNEng[jN];
                HCalHits += 1.;
            }
        }
    }


    bool ECalEpCut;
    bool HCalEpCut;

    if (trackMom.Mag() > 0.5) // Different cuts for low-momentum tracks
    {
        double cutValueEcal = eCalCutFunction->Eval(trackMom.Mag());
        double cutWidthEcal = eCalCutWidthFunction->Eval(trackMom.Mag());

        ECalEpCut = (ECalEnergy/trackMom.Mag() >= 0.0 && ECalEnergy/trackMom.Mag() <= cutValueEcal+cutWidthEcal);
    }
    else
    {
    ECalEpCut = false;
    }

    if (trackMom.Mag() > 0.5) // Different cuts for low-momentum tracks
    {
        double cutValueHcal = hCalCutFunction->Eval(trackMom.Mag());
        double cutWidthHcal = hCalCutWidthFunction->Eval(trackMom.Mag());

        HCalEpCut = (HCalEnergy/trackMom.Mag() >= cutValueHcal-cutWidthHcal && HCalEnergy/trackMom.Mag() <= cutValueHcal+cutWidthHcal);
    }
    else
    {
        HCalEpCut = false;
    }

    
    bool EcalHitCut = (ECalHits < 6);
    bool HcalHitCut = (HCalHits < 9);

    if (ECalEpCut && HCalEpCut && EcalHitCut && HcalHitCut)
    {
        return true;
    }
    else
    {
        return false;
    }


}