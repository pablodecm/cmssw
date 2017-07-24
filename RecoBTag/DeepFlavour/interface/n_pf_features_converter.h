#ifndef RecoSV_DeepFlavour_n_pf_features_converter_h
#define RecoSV_DeepFlavour_n_pf_features_converter_h

#include "RecoBTag/DeepFlavour/interface/deep_helpers.h"

namespace deep {

  template <typename CandidateType,
            typename JetType,
            typename NeutralCandidateFeaturesType>
  void n_pf_features_converter(const CandidateType * n_pf,
                               const JetType & jet,
                               const float & drminpfcandsv,
                               NeutralCandidateFeaturesType & n_pf_features) {

    n_pf_features.pt = n_pf->pt();
    n_pf_features.eta = n_pf->eta();
    n_pf_features.phi = n_pf->phi();
    n_pf_features.ptrel = deep::catch_infs_and_bound(n_pf->pt()/jet.pt(),
                                                     0,-1,0,-1);
    n_pf_features.erel = deep::catch_infs_and_bound(n_pf->energy()/jet.energy(),
                                                    0,-1,0,-1);
    n_pf_features.puppiw = n_pf->puppiWeight();
    n_pf_features.phirel = deep::catch_infs_and_bound(fabs(reco::deltaPhi(n_pf->phi(),jet.phi())),
                                                      0,-2,0,-0.5);
    n_pf_features.etarel = deep::catch_infs_and_bound(fabs(n_pf->eta()-jet.eta()),
                                                      0,-2,0,-0.5);
    n_pf_features.deltaR = deep::catch_infs_and_bound(reco::deltaR(*n_pf,jet),
                                                      0,-0.6,0,-0.6);
    n_pf_features.isGamma = 0;
    if(fabs(n_pf->pdgId())==22)  n_pf_features.isGamma = 1;
    n_pf_features.HadFrac = n_pf->hcalFraction();

    n_pf_features.drminsv = deep::catch_infs_and_bound(drminpfcandsv,
                                                       0,-0.4,0,-0.4);

  } 


}

#endif //RecoSV_DeepFlavour_n_pf_features_converter_h
