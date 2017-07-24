import FWCore.ParameterSet.Config as cms

from RecoBTag.DeepFlavour.DeepFlavourTagInfos_cfi import pfDeepFlavourTagInfos
from RecoBTag.DeepFlavour.DeepFlavourJetTags_cfi import pfDeepFlavourJetTags


pfDeepFlavourTaskNew = cms.Task(pfDeepFlavourTagInfos, pfDeepFlavourJetTags)
