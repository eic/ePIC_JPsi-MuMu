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

int ParentIndex(int simID, TTreeReaderArray<int>& partGenStat, TTreeReaderArray<int>& parentID, TTreeReaderArray<int>& partPdg, TTreeReaderArray<unsigned int>& partParb, TTreeReaderArray<unsigned int>& partPare)
{
    for (int k = partParb[simID]; k < partPare[simID]; k++)
    {
        int parentIndex = parentID[k];
        if ((partGenStat[parentIndex] == 1) && (partPdg[parentIndex] == 11 || TMath::Abs(partPdg[parentIndex]) == 13)) // Check if the parent is a electron or muon
        {
            return parentIndex;
        }
        else
        {
            for (int l = partParb[parentIndex]; l < partPare[parentIndex]; l++)
            {
                int grandparentIndex = parentID[l];
                if ((partGenStat[grandparentIndex] == 1) && (partPdg[grandparentIndex] == 11 || TMath::Abs(partPdg[grandparentIndex]) == 13)) // Check if the grandparent is a electron or muon
                {
                    return grandparentIndex;
                }
            }
        }
    }

    return -1; // Return -1 if no valid parent found

}