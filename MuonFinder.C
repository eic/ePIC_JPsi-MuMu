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

#include "DVMP_JPsi_Analysis.h"


bool IsMuon(TVector3 trackMom, int simuID,
            TTreeReaderArray<float>& EcalBarrelEng, TTreeReaderArray<float>& EcalEndcapPEng, TTreeReaderArray<float>& EcalEndcapNEng, 
            TTreeReaderArray<float>& HcalBarrelEng, TTreeReaderArray<float>& HcalEndcapPEng, TTreeReaderArray<float>& HcalEndcapNEng, 
            TTreeReaderArray<unsigned int>& simuAssocEcalBarrel, TTreeReaderArray<unsigned int>& simuAssocEcalEndcapP, TTreeReaderArray<unsigned int>& simuAssocEcalEndcapN,
            TTreeReaderArray<unsigned int>& simuAssocHcalBarrel, TTreeReaderArray<unsigned int>& simuAssocHcalEndcapP, TTreeReaderArray<unsigned int>& simuAssocHcalEndcapN)
{

    double ECalEnergy = 0.0;
    double HCalEnergy = 0.0;
    int ECalHits = 0;
    int HCalHits = 0;

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


    bool ECalEpCut = (ECalEnergy/trackMom.Mag() < 0.4 && ECalEnergy/trackMom.Mag() >= 0.);
    bool HCalEpCut = (HCalEnergy/trackMom.Mag() < 1 && HCalEnergy/trackMom.Mag() >= 0.05);
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