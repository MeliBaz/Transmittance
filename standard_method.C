void standard_method(TString input = "IR_Nikon-105_faces2_442nm_230926_1.root", int laserId = 4)
{
   vector<double> Mean_ratio_value;
   vector<double> Mean_ratio_ref;
   vector<double>Mean_value;
   vector<double>Mean_ref;
   vector<double>Double_ratio;
   vector<int> Position_x;
   vector<int> Position_y;
   vector<double> Ratio_value_gaus;
   vector<double> Ratio_ref_gaus;
   vector<double>Value;
   vector<double>Reference;

   
    struct Correction
    {
        string laser;
        double value;
        double err_sys;
        double err_stat;
    };

    Correction corr[] = {{"none", 1, 0, 0},
                         {"uv", 1, 0, 0},
                         {"lightuv", 1.03645846, 0.0, 0.0},
                         {"blue", 0.98172, 0.0, 0.00064},
                         {"lightblue", 0.974811, 0.0, 0.00064},
                         {"green", 0.959840, 0.0, 0.00084},
                         {"red", 0.95, 0.0, 0.00052},
    };
    

    Correction cor_mirror = corr[laserId];

    Double_t offsetMean_ref = +0.00073; // S1227-1010BR
    Double_t offsetMean_ref_err = 0.00034;

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

   TFile *file = new TFile(input, "READ");
    TTree *tree = (TTree *)file->Get("tree");
    Double_t reference, x_scan, y_scan, value;
    Int_t filter;
    tree->SetBranchAddress("reference", &reference);
    tree->SetBranchAddress("y_scan", &y_scan);
    tree->SetBranchAddress("x_scan", &x_scan);
    tree->SetBranchAddress("value", &value);  
    tree->SetBranchAddress("filter", &filter);
    Int_t numberOfEntries = tree->GetEntries();
    Double_t min_x = 1000;
    Double_t max_x = 0;
    Double_t min_y = 1000;
    Double_t max_y = 0;
    Double_t first_x, second_x = 0, first_y, second_y = 0;
    double_t step_x, step_y;
    Int_t nx, ny, p_filter;
    Double_t p_x, p_y;
    Double_t ratio_ref, ratio_value, double_ratio;
    Double_t mean_value,mean_ref;
    double mean_value_g, mean_reference_g;
    double v_max=0,v_min=5,first_v, second_v = 0,step_v,nv;
    double r_max=0,r_min=5,first_r, second_r = 0,step_r,nr;
                
    for (Int_t i = 0; i < numberOfEntries; i++)
    {
        tree->GetEntry(i);
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
    cout<<"x_min "<<min_x<< " "<<"y_min "<<min_y<<endl;
    cout<<"x_max "<<max_x<< " "<<"y_max "<<max_y<<endl;
   
    step_x = fabs(second_x - first_x);
    step_y = fabs(second_y - first_y);
    nx = ((max_x - min_x) / step_x)+1;
    ny = (max_y - min_y) / step_y +1;
    step_v = fabs(second_v - first_v);
    nv = ((v_max - v_min) / step_v) + 1;
    step_r = fabs(second_r - first_r);
    nr = ((r_max - r_min) / step_r) + 1;
 
    TH1F *hist_value = new TH1F("hist_value", "Value;Intensity;Entries", nv,v_min , v_max);
    TH1F *hist_reference = new TH1F("hist_reference", "Reference;Intensity;Entries", nr, r_min, r_max);
    TH2F *hist_ratio_map = new TH2F("hist_ratio_map", "Ratio_map;x [mm];y [mm]", nx, min_x-0.5*step_x, max_x+0.5*step_x, ny, min_y-0.5*step_y, max_y+0.5*step_y);
    TH1F *hist_x = new TH1F("Transmittance", "Transmittance;Intensity;Entries",150 , 0.9, 1.1);

   
     cout<<"nx "<<nx<<" "<<"ny "<<ny<<endl;

    for (Int_t i = 0; i < numberOfEntries; i++) 
    {
        tree->GetEntry(i);
        Position_x.push_back(x_scan);
        Position_y.push_back(y_scan); 
       
        //Bar in
            if (filter == 2) 
        {   
            ratio_value = (value - cor_offset.value) / (reference - offsetMean_ref);
            Mean_ratio_value.push_back(ratio_value);
            hist_value->Fill(value);
       
        }
        //Bar out
        if (filter == 0) 
        {
            ratio_ref = (value - cor_offset.value) / (reference - offsetMean_ref);
            Mean_ratio_ref.push_back(ratio_ref);
            hist_reference->Fill(reference);
         
        }
    }
  
   

    Int_t k=0;
    while(k<Mean_ratio_value.size() && k<Mean_ratio_value.size() ) 
    {
        double sum_value = 0;
        double sum_ref = 0;
        for (int j = k; j <k+10; j++) 
        {
            sum_value += Mean_ratio_value[j];
            sum_ref += Mean_ratio_ref[j];
            if (j==k+9){
                mean_value=sum_value/10.0;
                mean_ref=sum_ref / 10.0;
                Mean_value.push_back(mean_value);
                Mean_ref.push_back(mean_ref);
            }    
            
        }

      
        double_ratio=(mean_value/mean_ref)*cor_mirror.value;
        cout<<"Double_ratio: "<<double_ratio<<endl;
        cout<<"Position x : "<<Position_x[2*k+1]<<endl;
        cout<<"Position y : "<<Position_y[2*k+1]<<endl;
        hist_x->Fill(double_ratio);
        hist_ratio_map->Fill(Position_x[2*k+1], Position_y[2*k+1], double_ratio);
      

        k+=10;
    }  
    TCanvas *c5=new TCanvas();
    hist_x->Draw();

    TCanvas *c2=new TCanvas("c2"," ",400,200);
    c2->Divide(2,1);
    c2->cd(1);
    hist_ratio_map->Draw("colz text");
    c2->cd(2);
    hist_x->Draw();




}
