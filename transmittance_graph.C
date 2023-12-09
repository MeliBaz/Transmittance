double l=1200; //length of the quartz bar
double w=17; //width of the quartz bar
TF1 *refn = new TF1("refn","sqrt(1+0.6961663/(1-pow(0.0684043/(x*0.001),2))+0.4079426/(1-pow(0.1162414/(x*0.001),2))+0.8974794/(1-pow(9.896161/(x*0.001),2)))",150,700);
double reflections( double wavelength) {
    double n=refn->Eval(wavelength);
    double br_angle=TMath::ATan(n);
    double br_angle_inradiator=TMath::ASin((1/n)*sin(br_angle));
    double N=l*TMath::Tan(br_angle_inradiator)/w;
    cout<<wavelength<<" n "<<n<< " Alpha "<<br_angle_inradiator*TMath::RadToDeg()<< " br angle "<<br_angle*TMath::RadToDeg()<<" Reflections "<<N<<endl;
    return N;
}
double apsorbtion_length(double Wavelength){
    double L=500000.0/pow((442./Wavelength),4);
    return L;
}
void transmittance_graph(TString input="measurements/double_ratio_*.root"){
Double_t mean, mean_error,R;
Int_t Wavelength;
double L;
double N;
double delta_l;
double delta_R;

TChain ch("double_ratio");
ch.Add(input),
ch.SetBranchAddress("mean", &mean);
ch.SetBranchAddress("mean_error", &mean_error);
ch.SetBranchAddress("Wavelength", &Wavelength);
TGraphErrors* graph = new TGraphErrors();

for (Int_t i = 0; i < ch.GetEntries(); i++)
{   
    ch.GetEntry(i);
    N=reflections(Wavelength);
    L=apsorbtion_length(Wavelength*1.0);
    double br_angle=TMath::ATan(refn->Eval(Wavelength));
    R=TMath::Power(mean*exp(-l/L),1/N);
    //Calculating errors
    delta_l=sqrt(pow(1/cos(br_angle),2)+pow(1*sin(br_angle)*0.1/pow(cos(br_angle),2),2) );
    delta_R=sqrt(pow((1/N)*pow(mean,(1-N)/N)*exp(-l/(N*L))*mean_error,2)+pow(-pow(mean,1/N)*exp(-l/(L*N))*delta_l/(L*N),2) );
    cout<<"Mean : "<<mean<<" Mean error : "<<mean_error<<" Wavelength : "<<Wavelength<<" exp "<<exp(-l/L)<<" R : " <<R<<" delta l "<<delta_l<<" delta R "<<delta_R<<endl;
    graph->SetPoint(i,Wavelength,R);
    graph->SetPointError(i,0,delta_R);
    
}
// Plotting predictions from scalar scattering theory for surface roughness 5, 10, 15 angstroms
TF1 *func5 = new TF1("func5","1-pow(4* TMath::Pi()*0.5*cos(TMath::ATan(refn(x)))*refn(x)/x,2)",150,700);
TF1 *func10 = new TF1("func10","1-pow(4* TMath::Pi()*1*cos(TMath::ATan(refn(x)))*refn(x)/x,2)",150,700);
TF1 *func15 = new TF1("func15","1-pow(4* TMath::Pi()*1.5*cos(TMath::ATan(refn(x)))*refn(x)/x,2)",150,700);

TCanvas *c1=new TCanvas();
graph->SetTitle("IR_Nikon");
graph->SetMarkerColor(1);
graph->SetMarkerSize(1.5);
graph->SetMarkerStyle(21);
graph->Draw("ALP");
graph->GetYaxis()->SetRangeUser(0.995,1.0006);
graph->GetXaxis()->SetLimits(160,700);
graph->GetXaxis()->SetTitle("Wavelength (nm)");
graph->GetYaxis()->SetTitle("Reflectivity");
c1->SetGrid();

func15->SetLineColor(kRed);
func10->SetLineColor(kBlue);
func5->SetLineColor(kGreen);

TLegend *legend = new TLegend(0.7,0.1,0.9,0.3);
legend->SetTextSize(0.04);
legend->AddEntry(func5," 5 #AA","l");
legend->AddEntry(func10,"10 #AA","l");
legend->AddEntry(func15,"15 #AA","l");
legend->Draw("same");
func15->Draw("same");
func10->Draw("same");
func5->Draw("same");
c1->Print("TGraph.pdf");
c1->Print("TGraph.root");
}
