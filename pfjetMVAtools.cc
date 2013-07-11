#include <vector>
#include "$CMSSW_BASE/src/CMS2/NtupleMacros/CORE/CMS2.h"

vector <double> getGoodMVAs()
{
  vector <double> goodmvas;
  for( size_t ucjeti = 0; ucjeti < cms2.pfjets_p4().size(); ucjeti++) {   // uncorrected jets collection                                           
	for( size_t cjeti = 0; cjeti < cms2.pfjets_p4().size(); cjeti++) {   // corrected jets collection                                           
	  double deta = cms2.pfjets_p4().at(ucjeti).eta() - (cms2.pfjets_corL1FastL2L3().at(cjeti) * cms2.pfjets_p4().at(cjeti)).eta();
	  double dphi = acos(cos(cms2.pfjets_p4().at(ucjeti).phi() - (cms2.pfjets_corL1FastL2L3().at(cjeti) * cms2.pfjets_p4().at(cjeti)).phi()));
	  double dr = sqrt(deta*deta + dphi*dphi);
	  if (dr > 0.01) continue;
	  goodmvas.push_back(cms2.pfjets_full53xmvavalue().at(ucjeti));
	}
  }
  return goodmvas;
}
