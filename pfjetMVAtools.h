#include <vector>
#include "Math/VectorUtil.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

#ifndef PFJETMVATOOLS_H
#define PFJETMVATOOLS_H
bool sortByPFJetPt (const std::pair <LorentzVector, Int_t> &pfjet1, const std::pair<LorentzVector, Int_t> &pfjet2);
bool getGoodMVAs(vector <float> &goodmvas, std::string variable);
#endif
