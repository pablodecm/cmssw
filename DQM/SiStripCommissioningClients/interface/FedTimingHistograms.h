#ifndef DQM_SiStripCommissioningClients_FedTimingHistograms_H
#define DQM_SiStripCommissioningClients_FedTimingHistograms_H

#include "DQM/SiStripCommissioningClients/interface/CommissioningHistograms.h"
#include "DQM/SiStripCommissioningSummary/interface/FedTimingSummaryFactory.h"
#include "DQM/SiStripCommissioningAnalysis/interface/FedTimingAnalysis.h"

class MonitorUserInterface;

class FedTimingHistograms : public CommissioningHistograms {

 public:
  
  FedTimingHistograms( MonitorUserInterface* );
  virtual ~FedTimingHistograms();

  typedef SummaryHistogramFactory<FedTimingAnalysis> Factory;
  
  /** */
  void histoAnalysis( bool debug );

  /** */
  void createSummaryHisto( const sistrip::SummaryHisto&,
			   const sistrip::SummaryType&,
			   const std::string& top_level_dir,
			   const sistrip::Granularity& );

 protected:

  std::map<uint32_t,FedTimingAnalysis> data_;

  std::auto_ptr<Factory> factory_;
  
  const float optimumSamplingPoint_;
  float minDelay_;
  float maxDelay_; 
  uint32_t deviceWithMinDelay_;
  uint32_t deviceWithMaxDelay_;

};

#endif // DQM_SiStripCommissioningClients_FedTimingHistograms_H

