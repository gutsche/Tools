// system include files
#include <memory>
#include <cmath>
#include "TMath.h"
#include <algorithm>
#include <TLorentzVector.h>
#include "MetUtilities.h"

using namespace std;

MetUtilities::MetUtilities() { 
  
  //Tight Id
  fMVACut[0][0][0] =  0.5; fMVACut[0][0][1] = 0.6; fMVACut[0][0][2] = 0.6; fMVACut[0][0][3] = 0.9;
  fMVACut[0][1][0] = -0.2; fMVACut[0][1][1] = 0.2; fMVACut[0][1][2] = 0.2; fMVACut[0][1][3] = 0.6;
  fMVACut[0][2][0] =  0.3; fMVACut[0][2][1] = 0.4; fMVACut[0][2][2] = 0.7; fMVACut[0][2][3] = 0.8;
  fMVACut[0][3][0] =  0.5; fMVACut[0][3][1] = 0.4; fMVACut[0][3][2] = 0.8; fMVACut[0][3][3] = 0.9;
  //Medium id
  fMVACut[1][0][0] =  0.2; fMVACut[1][0][1] = 0.4; fMVACut[1][0][2] = 0.2; fMVACut[1][0][3] = 0.6;
  fMVACut[1][1][0] = -0.3; fMVACut[1][1][1] = 0. ; fMVACut[1][1][2] = 0. ; fMVACut[1][1][3] = 0.5;
  fMVACut[1][2][0] =  0.2; fMVACut[1][2][1] = 0.2; fMVACut[1][2][2] = 0.5; fMVACut[1][2][3] = 0.7;
  fMVACut[1][3][0] =  0.3; fMVACut[1][3][1] = 0.2; fMVACut[1][3][2] = 0.7; fMVACut[1][3][3] = 0.8;
  //Loose Id 
  fMVACut[2][0][0] =  0. ; fMVACut[2][0][1] =  0. ; fMVACut[2][0][2] =  0. ; fMVACut[2][0][3] = 0.2;
  fMVACut[2][1][0] = -0.4; fMVACut[2][1][1] = -0.4; fMVACut[2][1][2] = -0.4; fMVACut[2][1][3] = 0.4;
  fMVACut[2][2][0] =  0. ; fMVACut[2][2][1] =  0. ; fMVACut[2][2][2] =  0.2; fMVACut[2][2][3] = 0.6;
  fMVACut[2][3][0] =  0. ; fMVACut[2][3][1] =  0. ; fMVACut[2][3][2] =  0.6; fMVACut[2][3][3] = 0.2;

}

MetUtilities::~MetUtilities() {}

bool MetUtilities::passMVA(std::pair<LorentzVector,double> jet) { 
  int ptId = 0; 
  if(jet.first.pt() > 10 && jet.first.pt() < 20) ptId = 1;
  if(jet.first.pt() > 20 && jet.first.pt() < 30) ptId = 2;
  if(jet.first.pt() > 30                       ) ptId = 3;
  
  int etaId = 0;
  if(fabs(jet.first.eta()) > 2.5  && fabs(jet.first.eta()) < 2.75) etaId = 1; 
  if(fabs(jet.first.eta()) > 2.75 && fabs(jet.first.eta()) < 3.0 ) etaId = 2; 
  if(fabs(jet.first.eta()) > 3.0  && fabs(jet.first.eta()) < 5.0 ) etaId = 3; 
  if(jet.second  > fMVACut[2][ptId][etaId]) return true;
  return false;
}
MetUtilities::LorentzVector *MetUtilities::leadPt(std::vector<JetInfo> &iJets,bool iFirst) {
  double lPt0 = 0; double lPt1 = 0;
  LorentzVector *lLead=0; 
  LorentzVector *l2nd =0;
  for(int i0 = 0; i0 < int(iJets.size()); i0++) {
    if(iJets[i0].p4.pt() > lPt0) {lPt0 = iJets[i0].p4.pt(); lLead = &(iJets[i0].p4); continue;}
    if(iJets[i0].p4.pt() > lPt1) {lPt1 = iJets[i0].p4.pt(); l2nd  = &(iJets[i0].p4); continue;}
  }
  if(iFirst) return lLead;
  return l2nd;
}
int MetUtilities::NJets(std::vector<JetInfo> &iJets,double iPt) {
  int lNJet = 0; for(int i0 = 0; i0 < int(iJets.size()); i0++) if(iJets[i0].p4.pt() > iPt) lNJet++;
  return lNJet;
}
double MetUtilities::deltaR(LorentzVector &iVec1,LorentzVector &iVec2) { 
  double lDPhi = fabs(iVec1.phi()-iVec2.phi());
  if(lDPhi > 2.*TMath::Pi()-lDPhi) lDPhi =  2.*TMath::Pi()-lDPhi;
  double lDEta = iVec1.eta()-iVec2.eta();
  return sqrt(lDPhi*lDPhi+lDEta*lDEta);
}
void MetUtilities::cleanJets(std::vector<LorentzVector> &iVis,std::vector<JetInfo> &iJets) { 
  for(int i0   = 0; i0 < int(iJets.size()); i0++) { 
    for(int i1 = 0; i1 < int(iVis.size());  i1++) { 
      if(deltaR(iVis[i1],iJets[i0].p4) < 0.5) {
	iJets.erase (iJets.begin()+i0); i0--;
      }
    }
  }
}
std::pair<MetUtilities::LorentzVector,double> MetUtilities::TKMet(std::vector<std::pair<LorentzVector,double> > &iCands,double iDZ,int iLowDz) { 
  double            lSumEt = 0;
  LorentzVector     lVec;
  for(int i0 = 0; i0 < int(iCands.size()); i0++) {
    if(iCands[i0].second < 0   &&  iLowDz != 2) continue;
    if(iCands[i0].second > iDZ &&  iLowDz == 0) continue;
    if(iCands[i0].second < iDZ &&  iLowDz == 1) continue;
    lVec    -= iCands[i0].first;
    lSumEt  += iCands[i0].first.pt();
  }
  std::pair<LorentzVector,double> lTKMet(lVec,lSumEt);
  return lTKMet;
}
std::pair<MetUtilities::LorentzVector,double> MetUtilities::JetMet(std::vector<JetInfo> &iJets,bool iPassMVA) { 
  double            lSumEt = 0;
  LorentzVector           lVec;
  int lNPass = 0;
  for(int i0 = 0; i0 < int(iJets.size()); i0++) { 
    std::pair<LorentzVector,double> pMVAInfo(iJets[i0].p4,iJets[i0].mva);
    bool pPass =  passMVA(pMVAInfo);
    if( passMVA(pMVAInfo)  && !iPassMVA) continue;
    if(!passMVA(pMVAInfo)  &&  iPassMVA) continue;
    LorentzVector  pFullVec; pFullVec = iJets[i0].p4; //Full 4 vector
    //Now make the Neutral contributions
    TLorentzVector pTVec; pTVec.SetPtEtaPhiM(pFullVec.Pt()*iJets[i0].neutFrac,pFullVec.eta(),pFullVec.phi(),pFullVec.mass());
    LorentzVector  pVec ; pVec .SetCoordinates(pTVec.Px(),pTVec.Py(),pTVec.Pz(),pTVec.E());
    lVec    -= pVec;
    lSumEt  += pVec.Pt();
    lNPass++;
  }
  std::pair<LorentzVector,double> lJetMet(lVec,lSumEt);
  return lJetMet;
}
std::pair<MetUtilities::LorentzVector,double> MetUtilities::NoPUMet(std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<JetInfo> &iJets,double iDZ) { 
  double            lSumEt = 0;
  LorentzVector     lVec;
  std::pair<LorentzVector,double> lTKMet  = TKMet (iCands,iDZ,0);
  std::pair<LorentzVector,double> lJetMet = JetMet(iJets     ,true);
  lVec += lTKMet.first;     lSumEt += lTKMet .second;
  lVec += lJetMet.first;    lSumEt += lJetMet.second; 
  std::pair<LorentzVector,double> lNoPUMet(lVec,lSumEt);
  return lNoPUMet;
}
std::pair<MetUtilities::LorentzVector,double> MetUtilities::PUMet  (std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<JetInfo> &iJets,double iDZ) { 
  double            lSumEt = 0;
  LorentzVector     lVec;
  std::pair<LorentzVector,double> lTKMet  = TKMet (iCands,iDZ,   1);
  std::pair<LorentzVector,double> lJetMet = JetMet(iJets     ,false);
  lVec += lTKMet .first;    lSumEt += lTKMet .second;
  lVec += lJetMet.first;    lSumEt += lJetMet.second; 
  std::pair<LorentzVector,double> lPUMet(lVec,lSumEt);
  return lPUMet;
}
std::pair<MetUtilities::LorentzVector,double>  MetUtilities::PUCMet (std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<JetInfo> &iJets,double iDZ) { 
  double            lSumEt = 0;
  LorentzVector     lVec;
  std::pair<LorentzVector,double> lPFMet  = TKMet (iCands,iDZ,   2);
  std::pair<LorentzVector,double> lTKMet  = TKMet (iCands,iDZ,   1);
  std::pair<LorentzVector,double> lJetMet = JetMet(iJets     ,false);
  lVec += lPFMet .first;     lSumEt += lPFMet .second;
  lVec -= lTKMet .first;     lSumEt -= lTKMet .second; 
  lVec -= lJetMet.first;     lSumEt += lJetMet.second;  //The Sum Et line here was a bug in the training
  std::pair<LorentzVector,double> lPUCMet(lVec,lSumEt);
  return lPUCMet;
}
std::pair<MetUtilities::LorentzVector,double> MetUtilities::PFRecoil(double iSumEt,LorentzVector iVis,std::vector<std::pair<LorentzVector,double> > &iCands,double iDZ) { 
  std::pair<LorentzVector,double> lPFMet = TKMet(iCands,iDZ,2);
  LorentzVector lVec; lVec += lPFMet.first; lVec += iVis;
  double        lSumEt      = lPFMet.second-iSumEt; 
  std::pair<LorentzVector,double> lPFRecoil(lVec,lSumEt);
  return lPFRecoil;
}
std::pair<MetUtilities::LorentzVector,double> MetUtilities::TKRecoil(double iSumEt,LorentzVector iVis,std::vector<std::pair<LorentzVector,double> > &iCands,double iDZ) { 
  std::pair<LorentzVector,double> lTKMet = TKMet(iCands,iDZ,0);
  LorentzVector lVec; lVec += lTKMet.first; lVec += iVis;
  double lSumEt      = lTKMet.second-iSumEt; 
  std::pair<LorentzVector,double> lTKRecoil(lVec,lSumEt);
  return lTKRecoil;
}
std::pair<MetUtilities::LorentzVector,double> MetUtilities::NoPURecoil(double iSumEt,LorentzVector iVis,
								       std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<JetInfo> &iJets,
								       double iDZ) { 
  std::pair<MetUtilities::LorentzVector,double> lNPMet = NoPUMet(iCands,iJets,iDZ);
  LorentzVector lVec; lVec += lNPMet.first; lVec += iVis;
  double lSumEt      = lNPMet.second-iSumEt; 
  std::pair<LorentzVector,double> lNPRecoil(lVec,lSumEt);
  return lNPRecoil;
}
std::pair<MetUtilities::LorentzVector,double> MetUtilities::PUCRecoil(double iSumEt,LorentzVector iVis,
								      std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<JetInfo> &iJets,
								      double iDZ) { 
  std::pair<LorentzVector,double> lPUCMet = PUCMet(iCands,iJets,iDZ);
  LorentzVector lVec; lVec += lPUCMet.first; lVec += iVis;
  double lSumEt      = lPUCMet.second-iSumEt; 
  std::pair<LorentzVector,double> lPUCRecoil(lVec,lSumEt);
  return lPUCRecoil;
}
