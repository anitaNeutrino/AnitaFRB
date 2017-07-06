#include "Waterfall.h"
#include "AnitaDataset.h" 
#include "AnitaVersion.h" 
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

  std::vector<double> spectra[NUM_PHI][131];  //TODO; could theoretically upsample
  std::vector<double> heading; 
  std::vector<double> time; 

  //this will hold the resampled  spectra / times 
  TGraph* resampled[NUM_PHI][131]; 
  
  
  int start_pol = opt.use_hpol ? 0 : 1; 
  int end_pol = opt.use_vpol ? 1 : 0; 
  int npol = end_pol - start_pol+1; 

  FilterStrategy strat; 

  if (anita_version == 3) strat.addOperation(&alfa); 

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

      time.push_back(d.header()->triggerTime + d.header()->triggerTimeNs * 1e-9); 
      heading.push_back(d.gps()->heading); 

      size_t idx = heading.size(); 


      for (int phi = 0; phi < NUM_PHI; phi++) 
      {
        for (int  k = 0; k <131; k++) 
        {
          //only reserve an extra 5000 at a time to avoid blowing up too much in memory 
          if (spectra[phi][k].size() == spectra[phi][k].capacity())
          {
                spectra[phi][k].reserve(5000+spectra[phi][k].capacity()); 
          }
          spectra[phi][k].push_back(0); 
        }
      }


      //compute the spectra 
      
      FilteredAnitaEvent ev(d.useful(), &strat, d.gps(), d.header()); 

      //cutout likely paylaod blasts 
      double max_ratio_hpol, max_ratio_vpol; 
      ev.getMinMaxRatio(AnitaPol::kHorizontal, &max_ratio_hpol,0,0,0);
      if (max_ratio_hpol > opt.max_bottom_top_ratio) continue; 
      ev.getMinMaxRatio(AnitaPol::kVertical, &max_ratio_vpol,0,0,0); 
      if (max_ratio_vpol > opt.max_bottom_top_ratio) continue; 


      for (int phi = 0; phi < NUM_PHI; phi++)
      {
        for (int ring = 0; ring < 3; ring++)
        {
          for (int pol = start_pol; pol <= end_pol; pol++) 
          {
            const AnalysisWaveform * wf = ev.getFilteredGraph(ring * NUM_PHI + phi, AnitaPol::AnitaPol_t(pol)); 
            const TGraphAligned * pow = wf->power(); 
            for (int  k = 0; k <131; k++) 
            {
              spectra[phi][k][idx] +=  pow->GetY()[k]  / (3*npol); 
            }
          }
        }
      }
    }
  }

 
  //Now we have to resample each of our spectra to our desired time base 

  for (int phi = 0; phi < NUM_PHI; phi++)
  {
    for (int k = 0; k < 131; k++) 
    {
      TGraph g(time.size(), &time[0], &spectra[phi][k][0]); 
      resampled[phi][k] = (TGraph*) FFTtools::getInterpolatedGraphSparseInvert(&g, 1./opt.desired_Hz, (opt.end_time - opt.start_time+1) * opt.desired_Hz,  opt.waterfallInterpolationWidth); 

      //free the memory from this phi sectors raw values 
      std::vector<double>().swap(spectra[phi][k]);
    }
  }


  //phew, almost done 
  
  double ang_bin_width = 360 / opt.angular_bins; 
  TGraph gheading(time.size(), &time[0], &heading[0]); 
  gheading.SetBit(TGraph::kIsSortedX); 

  for (int ang = 0; ang < opt.angular_bins; ang++)
  {
    angles.push_back((ang+0.5)* ang_bin_width); 

    TString name; name.Form("waterfall_%d_%d_%d", opt.start_time, opt.end_time, (int) angles[ang]); 
    TString title; title.Form("Waterfall (bearing %g-%g)", angles[ang]-ang_bin_width/2, angles[ang] + ang_bin_width/2); 
    waterfalls.push_back(new TH2D(name.Data(), title.Data(), (opt.end_time - opt.start_time + 1) * opt.desired_Hz, opt.start_time, opt.end_time+1, 
                                                             131, 0, 1300)); 
    waterfalls[ang]->GetXaxis()->SetTimeDisplay(1); 
    waterfalls[ang]->GetXaxis()->SetTitle("time"); 
    waterfalls[ang]->GetYaxis()->SetTitle("f (MHz)"); 
    waterfalls[ang]->GetZaxis()->SetTitle("power"); 


    //take average of phi_sectors within width 
    //this could eventually do something smarter
   
    for (int ii = 1; ii < waterfalls[ang]->GetNbinsX(); ii++)
    {
      double head = gheading.Eval(waterfalls[ang]->GetXaxis()->GetBinCenter(ii)); 

      int nused = 0; 

      for (int phi  = 0; phi < NUM_PHI; phi++)
      {
        //TODO: check if the sign of the 45 degrees is right, otherwise everything will be off by 90 degrees 
        double center = FFTtools::wrap((phi) * 22.5 - 45-head,360); 

        if ( center >= angles[ang]-ang_bin_width/2 && center <= angles[ang] + ang_bin_width/2)
        {
          nused++; 
          for (int jj = 1; jj < waterfalls[ang]->GetNbinsY(); jj++)
          {
            waterfalls[ang]->SetBinContent(ii,jj, waterfalls[ang]->GetBinContent(ii,jj) + resampled[phi][jj-1]->GetY()[ii-1]); 
          }
        }
      }

      //normalize 
      for (int jj = 1; jj < waterfalls[ang]->GetNbinsY(); jj++)
      {
        waterfalls[ang]->SetBinContent(ii,jj, waterfalls[ang]->GetBinContent(ii,jj)/ nused); 
      }
    }
  }

  //delete resampled graphs we don't need anymore 
  for (int phi = 0; phi < NUM_PHI; phi++)
  {
    for (int k = 0; k < 131; k++) 
    {
      delete resampled[phi][k]; 
    }
  }


}


AnitaFRB::Waterfall::~Waterfall()
{

  for (size_t i = 0; i < nangles(); i++)
  {
    delete waterfalls[i]; 

  }

}
