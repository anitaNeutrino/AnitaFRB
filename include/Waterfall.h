#ifndef ANITA_FRB_WATERFALL_H
#define ANITA_FRB_WATERFALL_H

/* Makes appropriate waterfall style plots for FRB search. 
 *
 * This is a bit tricky and involves a few things: 
 *
 * 1) interpolating our event rate into bins of even sizes  (10 Hz or so?) 
 * 2) adjusting for heading / payload rotation 
 * 
 *
 */

#include "AnitaConventions.h" 
#include <vector>
#include "TH2.h" 



namespace AnitaFRB
{

  class WaterfallOptions 
  {
    public: 
    WaterfallOptions(unsigned start, unsigned end) 
      : start_time(start), end_time(end)
    {
      angular_bins = 8; 
      use_hpol = true; 
      use_vpol = true; 
      max_bottom_top_ratio = 2.5; 
      include_minbias = true; 
      include_rf = true; 
      min_freq_bin=20; 
      max_freq_bin=120;
      desired_Hz = 20; 
    }

    unsigned start_time; ///start time for this average (unix time, UTC). no default! 
    unsigned end_time; ///end time for this average (unix time, UTC). 
    int angular_bins; ///number of angular bins. This is in bearing (so basically assumes payload is stationary over the period we are interested in). You probably want this to be 4 8 or 16. 
    bool use_hpol; ///true to include hpol in average(default true) 
    bool use_vpol; ///true to use vpol in average (default true) 
    double max_bottom_top_ratio;  ///for payload blast rejection, default 1.5
    bool include_minbias; /// true (default) to include minbias triggers
    bool include_rf; /// true (default) to include rf triggers; 

    int min_freq_bin; //minimum frequency bin (multiple of 10 MHz) 
    int max_freq_bin; //maximum frequency bin (multiple of 10 MHz), not inclusive 

    int desired_Hz; //if non-zero, resamples to this
    /**
     * TODO: resampling options
     *
     */
  };

  class Waterfall : public TObject 
  {

    public: 


      Waterfall() { ; } 
      Waterfall (const WaterfallOptions & opts); 
      virtual ~Waterfall(); 


      size_t nangles() const { return angles.size(); } 
      double getAngle(int i) const { return angles[i]; }
      const TH2 * getSpectrogram(int i ) const { return waterfalls[i]; } 

    private: 

      std::vector<double> angles; 
      std::vector<TH2*> waterfalls; 

      ClassDef(Waterfall,1); 
  }; 



}

#endif
