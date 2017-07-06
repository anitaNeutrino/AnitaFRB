#include "Waterfall.h"
#include "AnitaDataset.h" 
#include "AnitaVersion.h" 
#include "TFile.h" 
#include "FFTtools.h"
#include "RFInterpolate.h" 
#include "Adu5Pat.h" 
#include "FilteredAnitaEvent.h" 
#include "RawAnitaHeader.h" 
#include "TH2.h" 
#include "BasicFilters.h"
#include "FilterStrategy.h" 

#include <stdio.h> 

ClassImp(AnitaFRB::Waterfall); 

typedef double waterfall_t; 
typedef TH2D waterfall_hist; 

const double dt = 1./2.6; 


static ALFAFilter alfa; 

AnitaFRB::Waterfall::Waterfall(const WaterfallOptions & opt) 
{

  //first step is to get the spectrum at each event. We will store each phi sector's average, and then appropriately partition
  // among the angular bins according to the heading. 
  if(opt.end_time < opt.start_time) 
  {
    fprintf(stderr,"End time is before start time. That's not so good is it? \n"); 
    return; 
  }

  if (!opt.use_hpol && !opt.use_vpol) 
  {

    fprintf(stderr,"You asked for neither hpol nor vpol. We don't have any other pols. \n"); 
    return; 
  }


  int anita_version = AnitaVersion::getVersionFromUnixTime(opt.start_time); 
  int start_run = AnitaDataset::getRunAtTime(opt.start_time); 
  int end_run = AnitaDataset::getRunAtTime(opt.end_time); 

  /* oh boy this potentially a LOT of memory! */ 

  int nfreq_bins = opt.max_freq_bin - opt.min_freq_bin; 

  /* This will hold the power in each angular bin prior to putting it in the histogram */ 
  std::vector<std::vector<std::vector<waterfall_t> > > spectra(opt.angular_bins, std::vector<std::vector<waterfall_t> >(nfreq_bins)); 

  std::vector<double> time; 

  /*Preallocate some memory will reuse a lot */ 
  std::vector<std::vector<waterfall_t> > working_power(NUM_PHI, std::vector<waterfall_t>(nfreq_bins)); 
  std::vector<waterfall_t> working_angle_power(nfreq_bins); 
  
  int start_pol = opt.use_hpol ? 0 : 1; 
  int end_pol = opt.use_vpol ? 1 : 0; 
  int npol = end_pol - start_pol+1; 

  FilterStrategy strat; 

  if (anita_version == 3) strat.addOperation(&alfa); 

  double ang_bin_width = 360 / opt.angular_bins; 
  for (int ang = 0; ang < opt.angular_bins; ang++)
  {
    angles.push_back((ang+0.5)* ang_bin_width); 
  }

 
  for (int run = start_run; run <= end_run; run++)
  {
    AnitaDataset d(run,false,WaveCalType::kDefault, AnitaDataset::ANITA_ROOT_DATA, AnitaDataset::kNoBlinding); 

    for (int i = 0; i < d.N(); i++) 
    {
      d.getEntry(i); 
      if (d.header()->triggerTime < opt.start_time) continue; //could be made more efficient, but doubtful that this will be a bottleneck
      if (d.header()->triggerTime > opt.end_time) break; 

      if (! opt.include_minbias && ( (d.header()->trigType & 1) == 0)) continue; 
      if (! opt.include_rf && ((d.header()->trigType & 1) == 1)) continue; 


     //compute the spectra 
      
      FilteredAnitaEvent ev(d.useful(), &strat, d.gps(), d.header()); 

      //cutout likely paylaod blasts 
      double max_ratio_hpol, max_ratio_vpol; 
      ev.getMinMaxRatio(AnitaPol::kHorizontal, &max_ratio_hpol,0,0,0);
      if (max_ratio_hpol > opt.max_bottom_top_ratio) continue; 
      ev.getMinMaxRatio(AnitaPol::kVertical, &max_ratio_vpol,0,0,0); 
      if (max_ratio_vpol > opt.max_bottom_top_ratio) continue; 

      time.push_back(d.header()->triggerTime + d.header()->triggerTimeNs * 1e-9); 

      double heading = d.gps()->heading; 

      if (time.size() % 100 == 0) printf("(read in %lu events)\n", time.size()); 
      


      for (int phi = 0; phi < NUM_PHI; phi++)
      {
        if (phi == 0) memset(&working_power[phi][0],0, sizeof(waterfall_t) * nfreq_bins); 

        for (int ring = 0; ring < 3; ring++)
        {
          for (int pol = start_pol; pol <= end_pol; pol++) 
          {
            const AnalysisWaveform * wf = ev.getFilteredGraph(ring * NUM_PHI + phi, AnitaPol::AnitaPol_t(pol)); 
            const TGraphAligned * pow = wf->power(); 
            for (int  k = opt.min_freq_bin; k <opt.max_freq_bin; k++) 
            {
              working_power[phi][k-opt.min_freq_bin] +=  pow->GetY()[k]  / (3.*npol); 
            }
          }
        }
      }


      for (int ang = 0; ang < opt.angular_bins; ang++) 
      {

        int nused = 0; 
        memset(&working_angle_power[0],0,sizeof(waterfall_t)*nfreq_bins); 

        /* figure out which phi sectors we need */ 
        for (int phi = 0; phi < NUM_PHI; phi++) 
        {
          double center = FFTtools::wrap((phi) * 22.5 - 45-heading,360); 
          if ( center >= angles[ang]-ang_bin_width/2 && center <= angles[ang] + ang_bin_width/2)
          {
            nused++; 
            for (int k = 0; k < nfreq_bins; k++)
            {
              working_angle_power[k] += working_power[phi][k]; 
            }

          }
       
        }


        /* store it*/ 
        for (int  k = 0; k <nfreq_bins; k++) 
        {
          //only reserve an extra 5000 at a time to avoid blowing up too much in memory 
          if (spectra[ang][k].size() == spectra[ang][k].capacity())
          {
                spectra[ang][k].reserve(5000+spectra[ang][k].capacity()); 
          }

          spectra[ang][k].push_back(working_angle_power[k]/nused); 
        }
      }
    }
  }

  //add a final time
  time.push_back(opt.end_time+1); 

 
  
  /* now make the histograms */ 

  for (int ang = 0; ang < opt.angular_bins; ang++)
  {

    TString name; name.Form("waterfall_%d_%d_%d", opt.start_time, opt.end_time, (int) angles[ang]); 
    TString title; title.Form("Waterfall (bearing %g-%g)", angles[ang]-ang_bin_width/2, angles[ang] + ang_bin_width/2); 

    int nbins = opt.desired_Hz ? (opt.end_time - opt.start_time+1) * opt.desired_Hz : time.size()-1; 

    waterfalls.push_back(new waterfall_hist(name.Data(), title.Data(), nbins, opt.start_time, opt.end_time+1, 
                                                             nfreq_bins, 10 * opt.min_freq_bin, 10 * opt.max_freq_bin)); 

    if (!opt.desired_Hz)
      waterfalls[ang]->GetXaxis()->Set(time.size()-1, &time[0]); 

    waterfalls[ang]->GetXaxis()->SetTimeDisplay(1); 
    waterfalls[ang]->GetXaxis()->SetTitle("time"); 
    waterfalls[ang]->GetYaxis()->SetTitle("f (MHz)"); 
    waterfalls[ang]->GetZaxis()->SetTitle("power (mW / antenna / event ?) "); 
    waterfalls[ang]->SetStats(0); 

   
    /*if don't resample, just  copy into hist */ 

    if (!opt.desired_Hz)
    {
      for (int jj = 1; jj <= waterfalls[ang]->GetNbinsY(); jj++)
      {
        for (int ii = 1; ii <= waterfalls[ang]->GetNbinsX(); ii++)
        {
          waterfalls[ang]->SetBinContent(ii,jj, spectra[ang][jj-1][ii-1]); 
        }
      }
    }

    /* otherwise, we have to find the right times */

    else
    {
      int idata = 0; 
      for (int ii = 1; ii <= waterfalls[ang]->GetNbinsX(); ii++)
      {
        int navg = 0; 
        while(time[idata] < waterfalls[ang]->GetXaxis()->GetBinLowEdge(ii+1) && idata < (int) time.size() )
        {
          for (int jj = 1; jj <= waterfalls[ang]->GetNbinsY(); jj++)
          {
            waterfalls[ang]->SetBinContent(ii,jj, waterfalls[ang]->GetBinContent(ii,jj) + spectra[ang][jj-1][idata]);  
          }
          idata++; 
          navg++; 
        }

        for (int jj = 1; jj <= waterfalls[ang]->GetNbinsY(); jj++)
        {
            waterfalls[ang]->SetBinContent(ii,jj, waterfalls[ang]->GetBinContent(ii,jj)/navg);  
        }
      }
    }

    //now free the memory from the spectrum we just filled
    std::vector<std::vector<waterfall_t> >().swap(spectra[ang]); 
  }

}


AnitaFRB::Waterfall::~Waterfall()
{

  for (size_t i = 0; i < nangles(); i++)
  {
    delete waterfalls[i]; 

  }

}
