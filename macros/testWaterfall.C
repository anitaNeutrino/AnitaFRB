void testWaterfall(int start_time = 1420070400, int length = 3600)
{

  TFile f("waterfall.root","RECREATE"); 
  AnitaFRB::WaterfallOptions opt(start_time, start_time + length); 
  AnitaFRB::Waterfall w(opt); 
  f.cd(); 
  w.Write("waterfall"); 
}
