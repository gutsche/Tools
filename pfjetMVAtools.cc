#include <vector>
#include "Math/VectorUtil.h"
#include "./CMS2.h"
// #include "$CMSSW_BASE/src/CMS2/NtupleMacros/CORE/CMS2.h"

#include "pfjetMVAtools.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std;
using namespace tas;

//usage
//in your looper, add these lines:
// vector <float> goodmvas;
// getGoodMVAs(goodmvas);

//note: the function getGoodMVAs returns true if it returns a different collection 
//of mva values. Otherwise, if there is no problem with the current set of mva values
//it will return false, and give you back the original set of mva values.

bool sortByPFJetPt (const LorentzVector &pfjet1, const LorentzVector & pfjet2)
{
  return pfjet1.pt() > pfjet2.pt();
}


//this function will return the index map of an unsorted jet collection to the sorted collection.
vector <Int_t> sortcorrectedjets (const vector <LorentzVector> & corrpfjets)
{
  
  vector <Int_t> indices;
  vector <LorentzVector> sortedjets = corrpfjets;
  Int_t tempind = 0;

  sort(sortedjets.begin(), sortedjets.end(), sortByPFJetPt);

  for( size_t sjeti = 0; sjeti < corrpfjets.size(); sjeti ++ ){//loop over sorted jets
	for( size_t cjeti = 0; cjeti < corrpfjets.size(); cjeti ++ ){//loop over unsorted jets
	  if( !(abs( sortedjets.at(sjeti).pt() - corrpfjets.at(cjeti).pt() ) > 0.1*(numeric_limits<float>::epsilon()) ) ){
		tempind = cjeti;
	  }
	}  
	indices.push_back(tempind);
  }
  return indices;
}

//this is the main function to fix the mva bug
bool getGoodMVAs(vector <float> &goodmvas, string variable)
{
  
  vector <float> mva_variable;
  if(       variable == "mvavalue"            ){ mva_variable = cms2.pfjets_mvavalue();
  }else if( variable == "full5xmvavalue"      ){ mva_variable = cms2.pfjets_full5xmvavalue();
  }else if( variable == "full53xmvavalue"     ){ mva_variable = cms2.pfjets_full53xmvavalue();
  }else if( variable == "full53xmva_nvtx"     ){ mva_variable = cms2.pfjets_full53xmva_nvtx();
  }else if( variable == "full53xmva_d0"       ){ mva_variable = cms2.pfjets_full53xmva_d0();
  }else if( variable == "full53xmva_dZ"       ){ mva_variable = cms2.pfjets_full53xmva_dZ();
  }else if( variable == "full53xmva_beta"     ){ mva_variable = cms2.pfjets_full53xmva_beta();
  }else if( variable == "full53xmva_betaStar" ){ mva_variable = cms2.pfjets_full53xmva_betaStar();
  }else if( variable == "full53xmva_nCharged" ){ mva_variable = cms2.pfjets_full53xmva_nCharged();
  }else if( variable == "full53xmva_dRMean"   ){ mva_variable = cms2.pfjets_full53xmva_dRMean();
  }else if( variable == "full53xmva_frac01"   ){ mva_variable = cms2.pfjets_full53xmva_frac01();
  }else if( variable == "full53xmva_frac02"   ){ mva_variable = cms2.pfjets_full53xmva_frac02();
  }else if( variable == "full53xmva_frac03"   ){ mva_variable = cms2.pfjets_full53xmva_frac03();
  }else if( variable == "full53xmva_frac04"   ){ mva_variable = cms2.pfjets_full53xmva_frac04();
  }else if( variable == "full53xmva_frac05"   ){ mva_variable = cms2.pfjets_full53xmva_frac05();
  }else{
	cout<<"variable not found. Check input. Exiting."<<endl;
	exit(99);
  }

  //if no bug is detected, returns the original collection of the mvas stored in the cms2 ntuple.
  if( cms2.pfjets_p4().size() == mva_variable.size() ) {
   
	goodmvas = mva_variable;
	return false;
   
  }else{
   
	vector <bool> isgoodindex;
	vector <LorentzVector> cjets;
	double deta = 0.0;
	double dphi = 0.0;
	double dr = 0.0;

	if( cms2.evt_isRealData() ){
	  for( size_t cjeti = 0; cjeti < cms2.pfjets_p4().size(); cjeti++) {   // corrected jets collection                                           
		cjets.push_back((double)cms2.pfjets_corL1FastL2L3residual().at(cjeti) * cms2.pfjets_p4().at(cjeti));
	  }
	}else{
	  for( size_t cjeti = 0; cjeti < cms2.pfjets_p4().size(); cjeti++) {   // corrected jets collection                                           
		cjets.push_back((double)cms2.pfjets_corL1FastL2L3().at(cjeti) * cms2.pfjets_p4().at(cjeti));
	  }
	}
	
	vector <Int_t> goodindices = sortcorrectedjets(cjets);
	for( size_t ucjeti = 0; ucjeti < cms2.pfjets_p4().size(); ucjeti++) {   // uncorrected jets collection                                           
	  for( size_t cjeti = 0; cjeti < cms2.pfjets_p4().size(); cjeti++) {   // corrected jets collection                                           
		
		//buggy method
		if( cms2.evt_isRealData() ){
		  if( abs( cms2.pfjets_area().at(ucjeti) - cms2.pfjets_area().at(goodindices.at(cjeti))) > numeric_limits<float>::epsilon() ) continue;
		  if( fabs( cms2.pfjets_p4().at(ucjeti).eta() - (cms2.pfjets_corL1FastL2L3residual().at(goodindices.at(cjeti)) * cms2.pfjets_p4().at(goodindices.at(cjeti))).eta()) > 0.01 ) continue;
		}else{
		  if( abs( cms2.pfjets_area().at(ucjeti) - cms2.pfjets_area().at(goodindices.at(cjeti))) > numeric_limits<float>::epsilon() ) continue;
		  if( fabs( cms2.pfjets_p4().at(ucjeti).eta() - (cms2.pfjets_corL1FastL2L3().at(goodindices.at(cjeti)) * cms2.pfjets_p4().at(goodindices.at(cjeti))).eta()) > 0.01 ) continue;
		}
		
		//fix
		if( cms2.evt_isRealData() ){
		  deta = cms2.pfjets_p4().at(ucjeti).eta() - (cms2.pfjets_corL1FastL2L3residual().at(goodindices.at(cjeti)) * cms2.pfjets_p4().at(goodindices.at(cjeti))).eta();
		  dphi = acos(cos(cms2.pfjets_p4().at(ucjeti).phi() - (cms2.pfjets_corL1FastL2L3residual().at(goodindices.at(cjeti)) * cms2.pfjets_p4().at(goodindices.at(cjeti))).phi()));
		  dr = sqrt(deta*deta + dphi*dphi);
		}else{
		  deta = cms2.pfjets_p4().at(ucjeti).eta() - (cms2.pfjets_corL1FastL2L3().at(goodindices.at(cjeti)) * cms2.pfjets_p4().at(goodindices.at(cjeti))).eta();
		  dphi = acos(cos(cms2.pfjets_p4().at(ucjeti).phi() - (cms2.pfjets_corL1FastL2L3().at(goodindices.at(cjeti)) * cms2.pfjets_p4().at(goodindices.at(cjeti))).phi()));
		  dr = sqrt(deta*deta + dphi*dphi);
		}
		
		if (dr > 0.01){
		  isgoodindex.push_back(false);
		}else{
		  isgoodindex.push_back(true);
		}
	  }
	}

	//fill the new mva values
	for( size_t mvai = 0; mvai < mva_variable.size(); mvai++ ){
	  if( isgoodindex.at(mvai) ) goodmvas.push_back(mva_variable.at(mvai));	
	}
  
	//still possible that the fix picks up less events than the fix in cmssw
	//This behavior was not seen by me, but just in case this line here will 
	// prevent the code from crashing and return the original mva collection.
	if( goodmvas.size() >= cms2.pfjets_p4().size() ){
	  return true;  
	}else{
	  cout<<"new mva values vector size smaller than pfjets collection size."<<endl;
	  cout<<"returning old mva collection."<<endl;
	  goodmvas.clear();
	  goodmvas = mva_variable;
	  return false;
	}
  }
}
