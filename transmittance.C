void transmittance(TString input=" ",int laserId=-1){ 
  if (laserId<0){
    
    if (input.Contains("266nm"))  laserId=1;
    if (input.Contains("325nm"))  laserId=2;
    if (input.Contains("405nm"))  laserId=3;
    if (input.Contains("442nm"))  laserId=4;
    if (input.Contains("532nm"))  laserId=5;
    if (input.Contains("635nm"))  laserId=6;
  }
  int Wavelength[] = {0, 266, 325, 405, 442, 532, 635};
  double mean_value[6];
  struct Correction {
    string laser;
    double value;
    double err_sys;
    double err_stat;
  };
  //Mirror corrections for different wavelengths
  Correction corr[] = {
    {"none", 1, 0, 0},
    {"uv", 0.952951724, 0, 0},
    {"lightuv", 1.03645846, 0.0, 0.0},
    {"blue", 0.98172, 0.00064, 0},
    {"lightblue",0.97574, 0.01, 0},
    {"green", 0.959840,0.00084 ,0 },
    {"red", 0.95, 0.00052,0 },
  };
  //Mirror corrections from Marvin
    /*Correction corr[] = {
    {"none", 1, 0, 0},
    {"uv", 1.02313, 0, 0},
    {"lightuv", 1.01354, 0.0, 0.0},
    {"blue", 0.99835, 0.0, 0.00064},
    {"lightblue", 0.98173, 0.0, 0.00064},
    {"green", 0.96878, 0.0, 0.00084},
    {"red", 0.95227, 0.0, 0.00052},
  };*/
  Correction cor_mirror = corr[laserId];

  //Reference diode offset correction  
  double offsetMean_ref = +0.00073; // S1227-1010BR
  double offsetMean_ref_err = 0.00034;

  //Measurement diode offset correction  
  Correction corr_offset[] = {
    {"none", 0, 0, 0},
    {"uv", -0.00014, 0.00034, 0},
    {"lightuv", 0, 0.0, 0.0},
    {"blue", -0.00016, 0.00034, 0},
    {"lightblue", 0, 0.0, 0.0},
    {"green", -0.00009, 0.00037, 0},
    {"red", -0.00009, 0.00037, 0},
  };
  Correction cor_offset = corr_offset[laserId];

//Read data from the root file
  TFile *file = new TFile(input, "READ");
  TTree *tree = (TTree *)file->Get("tree");
  double reference, x_scan, y_scan, value, mreference, mvalue;
  int filter;
  tree->SetBranchAddress("reference", &reference);
  tree->SetBranchAddress("y_scan", &y_scan);
  tree->SetBranchAddress("x_scan", &x_scan);
  tree->SetBranchAddress("value", &value);
  tree->SetBranchAddress("filter", &filter);

int numberOfEntries = tree->GetEntries();

  double min_x = 1000, max_x = 0, min_y = 1000, max_y = 0;
  double first_x, second_x = 0, first_y, second_y = 0, step_x, step_y;
  int nx, ny, p_filter;
  double p_x, p_y;
  double ratio_ref, ratio_value, double_ratio;
  double v_max=0,v_min=5,first_v, second_v = 0,step_v,nv;
  double r_max=0,r_min=5,first_r, second_r = 0,step_r,nr;
 
 for (int i = 0; i < numberOfEntries; i++) {
    tree->GetEntry(i);
  
    //Find the range of data
    if (x_scan < min_x) min_x = x_scan;
    if (x_scan > max_x) max_x = x_scan;
    if (y_scan < min_y) min_y = y_scan;
    if (y_scan > max_y) max_y = y_scan;
    if (i == 0) first_x = x_scan;
    if (second_x == 0 && first_x != x_scan) second_x = x_scan;
    if (i == 0) first_y = y_scan;
    if (second_y == 0 && first_y != y_scan) second_y = y_scan;
    if (value < v_min) v_min = value;
    if (value> v_max) v_max= value;
    if (i == 0) first_v = value;
    if (second_v == 0 && first_v != value) second_v = value;
    if (reference < r_min) r_min = reference;
    if (reference> r_max) r_max= reference;
    if (i == 0) first_r = reference;
    if (second_r == 0 && first_r != reference) second_r = reference;

  }
  cout << "x_min " << min_x << " " << "y_min " << min_y << endl;
  cout << "x_max " << max_x << " " << "y_max " << max_y << endl;

//Calculate the number of bins
  step_x = fabs(second_x - first_x);
  step_y = fabs(second_y - first_y);
  nx = ((max_x - min_x) / step_x) + 1;
  ny = (max_y - min_y) / step_y + 1;
  step_v = fabs(second_v - first_v);
  nv = ((v_max - v_min) / step_v) + 1;
  step_r = fabs(second_r - first_r);
  nr = ((r_max - r_min) / step_r) + 1;
  
  TH1F *hist_value = new TH1F("hist_value", "Value;Intensity;Entries", nv,v_min , v_max);
  TH1F *hist_reference = new TH1F("hist_reference", "Reference;Intensity;Entries", nr, r_min, r_max);
  TH1F *hist_ratio_ref = new TH1F("hist_ratio_ref", "Ratio_Ref;Intensity;Entries", 150, 0.8, 2);
  TH1F *hist_ratio_value = new TH1F("hist_ratio_value", "Ratio_Value;Intensity;Entries", 100, 0.8, 2);
  TH1F *hist_double_ratio = new TH1F("Transmittance", "Transmittance;Intensity;Entries", 150, 0.8, 1.1);
  TH2F *hist_ratio_map = new TH2F("Ratio map", "Ratio map;x [mm];y [mm]", nx, min_x - 0.5 * step_x,max_x + 0.5 * step_x, ny, min_y - 0.5 * step_y, max_y + 0.5 * step_y);

  for (int i = 240; i < numberOfEntries; i++) {
    tree->GetEntry(i);

    //Define the initial positions and the filter inital value
    if (i == 0) {
      p_x = x_scan;
      p_y = y_scan;
      p_filter = filter;
    }
    //Calculate the double ratio for every 10 measurements
    if (x_scan != p_x || y_scan != p_y || filter != p_filter || i == numberOfEntries - 1) {
      //Fit histogram value with the Gaussian function and obtain mean value
      auto r = hist_value->Fit("gaus", "SQL");
      mvalue = r->Parameter(1);
      //mvalue=hist_value->GetMean();
      hist_value->Draw();
      //gPad->Update();
      //gPad->WaitPrimitive();

      //Fit histogram reference with the Gaussian function and obtain mean value
      r = hist_reference->Fit("gaus", "SQL");
      mreference = r->Parameter(1);
      //mreference=hist_reference->GetMean();
      //gPad->Update();
      //gPad->WaitPrimitive();

      //Calulate the ratio of value to reference based on the bar in (filter=2) and bar out (filter=0) filters
      if(p_filter == 0) {
        ratio_ref = (mvalue - cor_offset.value) / (mreference - offsetMean_ref);
        hist_ratio_ref->Fill(ratio_ref);
      }

      if(p_filter == 2) {
        ratio_value = (mvalue - cor_offset.value) / (mreference - offsetMean_ref);
        double_ratio = (ratio_value / ratio_ref) * (cor_mirror.value);	
  
        cout << p_filter<<" mvalue " << mvalue<<" "<<mreference << " double_ratio " << double_ratio << endl;
  
        hist_double_ratio->Fill(double_ratio);
        hist_ratio_map->Fill(p_x, p_y, double_ratio);
        hist_ratio_value->Fill(ratio_value);
    
      }
      //Reset histogram value and reference for every 10 measurements
      hist_value->Reset();
      hist_reference->Reset();

      p_x = x_scan;
      p_y = y_scan;
      p_filter = filter;
    }

    hist_reference->Fill(reference);
    hist_value->Fill(value);
  }

  //Fit histogram for double ratio with the Gaussian function and obtain mean value and mean error
  gStyle->SetOptFit(1);
  auto r = hist_double_ratio->Fit("gaus", "SL");
  double mean = r->Parameter(1);
  double mean_error = r->ParError(1);
  mean_error=sqrt(pow(mean_error,2)+pow(cor_mirror.err_sys,2));
  //cout<<"Corr "<<1/mean<<endl;

  //Create a file with data from the double ratio histogram to plot TGraphErrors (Double_ratio_graph.C) 
  TFile *output_data = new TFile(Form("double_ratio_%d.root", laserId), "recreate");
  TTree *otree = new TTree("double_ratio", "double_ratio");
  otree->Branch("mean", &mean);
  otree->Branch("mean_error", &mean_error);
  otree->Branch("Wavelength", &Wavelength[laserId]);
  otree->Fill();
  output_data->Write();

 
  hist_ratio_map->GetXaxis()->SetRangeUser(65, 73);
  hist_ratio_map->SetMinimum(0.5);
  TCanvas *c3 = new TCanvas("c3", "Raio reference", 800, 600);
  hist_ratio_ref->Draw();
  TCanvas *c4 = new TCanvas("c4", "Ratio value", 800, 600);
  hist_ratio_value->Draw();
  TCanvas *c5 = new TCanvas("c5", "Double ratio without fit", 800, 600);
  hist_double_ratio->Draw();
  TCanvas *c6 = new TCanvas("c6", "Ratio map", 800, 600);
  hist_ratio_map->Draw("colz text");
 
  TString folder = "hists/";
  gSystem->Exec("mkdir -p " +folder);
  c6->Print(folder + Form("double_ratio_%d.png", laserId));
  c3->Print(folder + Form("hist_ratio_ref_%d.root", laserId));
  c4->Print(folder + Form("hist_ratio_value_%d.root", laserId));
  c5->Print(folder + Form("hist_double_ratio_without_fit%d.root", laserId));
  c6->Print(folder + Form("hist_ratio_map_%d.root", laserId));

  c3->Print(folder + Form("hist_ratio_ref_%d.png", laserId));
  c4->Print(folder + Form("hist_ratio_value_%d.png", laserId));
  c5->Print(folder + Form("hist_double_ratio_%d.png", laserId));
  c6->Print(folder + Form("hist_ratio_map_%d.png", laserId));

  TCanvas *c11 = new TCanvas("c11", input, 400, 200);
  c11->SetTitle("Method with Gaussian fit");
  c11->Divide(2, 1);
  c11->cd(1);
  hist_ratio_map->Draw("colz text");
  c11->cd(2);
  hist_double_ratio->Draw();

  c11->Print("transmittance_with_fit.pdf[");
  c11->Print("transmittance_with_fit.pdf");
  c11->Print("transmittance_with_fit.pdf]");
  c3->Print("hists_transmittance.pdf[");
  c3->Print("hists_transmittance.pdf");
  c4->Print("hists_transmittance.pdf");
  c4->Print("hists_transmittance.pdf]");


}
