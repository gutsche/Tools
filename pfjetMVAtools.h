#include <vector>
#include "Math/VectorUtil.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

#ifndef PFJETMVATOOLS_H
#define PFJETMVATOOLS_H
bool sortByPFJetPt (const LorentzVector &pfjet1, const LorentzVector & pfjet2);
std::vector <Int_t> sortcorrectedjets (const std::vector <LorentzVector> & corrpfjets);
bool getGoodMVAs(std::vector <float> &goodmvas);
#endif
