#pragma once

// Define 3 and 4-momentum vectors for particles

inline TVector3 recoTrackMom[3];
inline ROOT::Math::PxPyPzEVector recoTrack4Mom[3];
inline int recoTrackIndex[3] = {0,0,0}; // Indices to keep track of which reconstructed track corresponds to which particle (e, mu+, mu-)
inline int recoTrackTruePID[3] = {0,0,0}; // True PDG codes for the reconstructed tracks

inline TVector3 beamEMom;
inline ROOT::Math::PxPyPzEVector beamE4Mom;
inline TVector3 scatEMomT;
inline ROOT::Math::PxPyPzEVector scatE4MomT;
inline TVector3 scatEMomR;
inline ROOT::Math::PxPyPzEVector scatE4MomR;

inline TVector3 beampMom;
inline ROOT::Math::PxPyPzEVector beamp4Mom;
inline TVector3 scatpMomT;
inline ROOT::Math::PxPyPzEVector scatp4MomT;
inline TVector3 scatpMomR;
inline ROOT::Math::PxPyPzEVector scatp4MomR;

inline TVector3 muPlusMomT;
inline ROOT::Math::PxPyPzEVector muPlus4MomT;
inline TVector3 muPlusMomR;
inline ROOT::Math::PxPyPzEVector muPlus4MomR;
inline TVector3 muMinusMomT;
inline ROOT::Math::PxPyPzEVector muMinus4MomT;
inline TVector3 muMinusMomR;
inline ROOT::Math::PxPyPzEVector muMinus4MomR;

inline ROOT::Math::PxPyPzEVector muPlus4MomR_Boosted;
inline ROOT::Math::PxPyPzEVector muMinus4MomR_Boosted;

inline TVector3 JPsiMomT;
inline ROOT::Math::PxPyPzEVector JPsi4MomT;
inline TVector3 JPsiMomR;
inline ROOT::Math::PxPyPzEVector JPsi4MomR;

inline double Q2_truth, Q2_DA, Q2_JB, Q2_e, Q2_sigma;
inline double t_truth, t_eXBABE, t_eXPT, t_eX, t_BABE;
inline double x_truth, x_DA, x_JB, x_e, x_sigma;
inline double y_truth, y_DA, y_JB, y_e, y_sigma;

inline double crossingAngle = 0.025; // 25 mrad crossing angle

inline double pMass = 0.93827208816; 
inline double eMass = 0.000510998950; 
inline double muMass = 0.1056583745;