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


bool IsElectron(TVector3 trackMom, int simuID,
            TTreeReaderArray<float>& EcalBarrelEng, TTreeReaderArray<float>& EcalEndcapPEng, TTreeReaderArray<float>& EcalEndcapNEng, 
            TTreeReaderArray<unsigned int>& simuAssocEcalBarrel, TTreeReaderArray<unsigned int>& simuAssocEcalEndcapP, TTreeReaderArray<unsigned int>& simuAssocEcalEndcapN)
{

    double ECalEnergy = 0.0;

    if (EcalBarrelEng.GetSize() == simuAssocEcalBarrel.GetSize())
    {
        for (int jB = 0; jB < simuAssocEcalBarrel.GetSize(); jB++) // Look for associations in the Ecal Barrel
        {
            if (simuAssocEcalBarrel[jB] == simuID)
            {
                ECalEnergy += EcalBarrelEng[jB];
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
            }
        }
    }



    bool CalEpCut = (ECalEnergy/trackMom.Mag() < 1.2 && ECalEnergy/trackMom.Mag() > 0.9);

    if (CalEpCut)
    {
        return true;
    }
    else
    {
        return false;
    }


}