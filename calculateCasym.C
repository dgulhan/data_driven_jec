#include "relativeResponse_MC.C"

void calculateCasym(){
TH1D::SetDefaultSumw2();
double etacut=3;
double phicut=2.5;
int docorrection =1;
char *mode="pPb";//"pPb","pp","Pbp";
char *algo="akPu3PF";
TString infile;
char *dataset="jet40";//jet80,MB;
char *binning="double_hcalbins";//"hcalbins"
double alphacut_low=0;
double alphacut_high=0.2;
int ptlow = 40;
int pthigh = 300;

TFile *f_mc = new TFile(Form("ntuples/ntuple_relativeResponse_MC_%s_eta%d_algo%s.root",mode,(int)etacut,algo));
TFile *f_data = new TFile(Form("ntuples/ntuple_relativeResponse_%s_eta%d_corrected%d_%s_%s.root",mode,(int)etacut,docorrection,dataset,algo));

TTree *nt_mc = (TTree*)f_mc->Get("ntdijet");
TTree *nt_data = (TTree*)f_data->Get("ntdijet");
TTree *nt_corr;
if(docorrection) nt_corr= (TTree*)f_data->Get("ntdijet_corr");
cout<<5<<endl;

int netabins;
if(binning=="double_hcalbins") netabins=28;
if(binning=="hcalbins") netabins=58;
double etabins[netabins];
double etabins_hcalbins[]= {-3, -2.853,
                     -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830,
                     -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
                     -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609,
                     -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,  0.000,
                      0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,
                      0.696,  0.783,  0.879,  0.957,  1.044,  1.131,  1.218,
                      1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,
                      1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  2.853,
                      3};
double etabins_double_hcalbins[]= {-3,-2.500,  -2.172, -1.930,-1.740,  -1.566, -1.392,  -1.218,-1.044,  -0.879,  -0.696,-0.522,  -0.348,  -0.174,  0.000,0.174,  0.348,    0.522,0.696,   0.879,   1.044,   1.218,1.392,   1.566,   1.740,1.930,   2.172,   2.500,   3};

                      
for(int i=0;i<netabins+1;i++){
 if(binning=="hcalbins"){
  cout<<1<<endl;
  etabins[i]=etabins_hcalbins[i];
 }
 if(binning=="double_hcalbins"){
  etabins[i]=etabins_double_hcalbins[i];
  cout<<etabins[i]<<endl;
 }
}

TProfile *B_mc = new TProfile("B_mc","",netabins,etabins);
  cout<<3<<endl;

TProfile *B_data = new TProfile("B_data","",netabins,etabins);
  cout<<4<<endl;

TProfile *B_corr = new TProfile("B_corr","",netabins,etabins);
  cout<<5<<endl;

nt_mc->Draw("(ptProbe-ptBarrel)/ptAverage:etaProbe>>B_mc",Form("weight*(ptAverage>%d && ptAverage<%d && alpha>=%.2f && alpha<%.2f && dphi>2.5)",ptlow,pthigh,alphacut_low,alphacut_high));
  cout<<6<<endl;

nt_data->Draw("(ptProbe-ptBarrel)/ptAverage:etaProbe>>B_data",Form("ptAverage>%d && ptAverage<%d && alpha>=%.2f && alpha<%.2f  && dphi>2.5",ptlow,pthigh,alphacut_low,alphacut_high));
  cout<<7<<endl;

TH1D * etaprobe_mc = new TH1D("etaprobe_mc","",netabins,etabins);
TH1D * etaprobe_data = new TH1D("etaprobe_data","",netabins,etabins);
TH1D * etaprobe_corr = new TH1D("etaprobe_corr","",netabins,etabins);
nt_mc->Draw("etaProbe>>etaprobe_mc",Form("weight*(ptAverage>%d && ptAverage<%d && alpha>=%.2f && alpha<%.2f && dphi>2.5)",ptlow,pthigh,alphacut_low,alphacut_high));
nt_data->Draw("etaProbe>>etaprobe_data",Form("ptAverage>%d && ptAverage<%d && alpha>=%.2f && alpha<%.2f && dphi>2.5",ptlow,pthigh,alphacut_low,alphacut_high));

TH1D * R_mc = new TH1D("R_mc","",netabins,etabins);
TH1D * R_data = new TH1D("R_data","",netabins,etabins);
TH1D * R_corr;
if(docorrection){
 nt_corr->Draw("(ptProbe-ptBarrel)/ptAverage:etaProbe>>B_corr",Form("ptAverage>%d && ptAverage<%d && alpha>=%.2f && alpha<%.2f && dphi>2.5",ptlow,pthigh),"prof");
 nt_corr->Draw("etaProbe>>etaprobe_corr",Form("ptAverage>%d && ptAverage<%d && alpha>=%.2f && alpha<%.2f && dphi>2.5",ptlow,pthigh));
 R_corr = new TH1D("R_corr","",netabins,etabins);
}

for(int i=0; i<netabins; i++){
// for(int i=0; i<100; i++){
 double Bval_mc = B_mc->GetBinContent(i+1);
 double Bval_data = B_data->GetBinContent(i+1);
 double Rerr_mc = 4*B_mc->GetBinError(i+1)/pow((2-Bval_mc),2);
 double Rerr_data = 4*B_data->GetBinError(i+1)/pow((2-Bval_data),2);
 R_mc->SetBinContent(i+1,(2+Bval_mc)/(2-Bval_mc));
 R_data->SetBinContent(i+1,(2+Bval_data)/(2-Bval_data));
 R_mc->SetBinError(i+1,Rerr_mc);
 R_data->SetBinError(i+1,Rerr_data);
 R_mc->SetMarkerColor(2);
 R_mc->SetLineColor(2);

 if(docorrection){
  double Bval_corr = B_corr->GetBinContent(i+1);
  double Rerr_corr = 4*B_corr->GetBinError(i+1)/pow((2-Bval_corr),2);
  R_corr->SetBinContent(i+1,(2+Bval_corr)/(2-Bval_corr));
  R_corr->SetBinError(i+1,Rerr_corr);
  R_corr->SetMarkerColor(kBlue);
  R_corr->SetLineColor(kBlue);
  R_corr->SetMarkerStyle(24);
 }
}
R_data->GetXaxis()->SetTitle("#eta_{probe}");
R_data->GetXaxis()->CenterTitle();
R_data->GetYaxis()->SetTitle("Relative Response");
R_data->GetYaxis()->CenterTitle();


TH1D *C_asym = (TH1D*)R_mc->Clone("C_asym");
C_asym->Divide(R_data);
C_asym->SetMaximum(1.25);
C_asym->SetMinimum(0.9);
C_asym->GetXaxis()->SetTitle("#eta_{probe}");
C_asym->GetXaxis()->CenterTitle();
C_asym->GetYaxis()->SetTitle("R_{MC}/R_{data}");
C_asym->GetYaxis()->CenterTitle();

TLegend *t3=new TLegend(0.55,0.8,0.9,0.93);
t3->SetFillColor(0);
t3->SetBorderSize(0);
t3->SetFillStyle(0);
t3->SetTextFont(63);
t3->SetTextSize(18); 
t3->AddEntry(R_data,"uncorrected","p");
if(docorrection) t3->AddEntry(R_corr,"corrected","p");
t3->AddEntry(R_mc,"PYTHIA + HIJING","p");

TCanvas *c1 = new TCanvas("c1","",600,600);
R_data->SetMaximum(1.25);
R_data->SetMinimum(0.75);
R_data->Draw();
R_mc->Draw("Same");
if(docorrection) R_corr->Draw("Same");
t3->Draw("Same");
drawText("anti-k_{T} PF R=0.3",0.2,0.90);
drawText("|#eta_{barrel}|<1.3",0.2,0.84);
drawText("|#Delta#phi_{1,2}|>2.5",0.2,0.77);
drawText(Form("%d<p_{T,ave}<%d GeV/c",ptlow,pthigh),0.2,0.71);
c1->SaveAs(Form("plots/R_%s_%s_%s_%s.png",mode,binning,algo,dataset));
c1->SaveAs(Form("plots/R_%s_%s_%s_%s.pdf",mode,binning,algo,dataset));

TCanvas *c2 = new TCanvas("c2","",600,600);
C_asym->Draw();
c2->SaveAs(Form("plots/Casym_%s_%s_%s_%s.png",mode,binning,algo,dataset));

TH1D * R_data_MC =(TH1D*)R_data->Clone("R_data_MC");
R_data_MC->Divide(R_mc);
TH1D * R_data_MC_corr;

if(docorrection){
 R_data_MC_corr=(TH1D*)R_corr->Clone("R_data_MC_corr");
 R_data_MC_corr->Divide(R_mc); 
}

TLegend *t4=new TLegend(0.55,0.8,0.9,0.93);
t4->SetFillColor(0);
t4->SetBorderSize(0);
t4->SetFillStyle(0);
t4->SetTextFont(63);
t4->SetTextSize(18); 
t4->AddEntry(R_data_MC,"uncorrected","p");
if(docorrection) t4->AddEntry(R_data_MC_corr,"corrected","p");

TLine l(-3,1,3,1);
TCanvas *c3 = new TCanvas("c3","",600,600);
 
R_data_MC->GetYaxis()->SetTitle("Data/MC");
R_data_MC->Draw();
if(docorrection) R_data_MC_corr->Draw("Same");
l.Draw("Same");
t4->Draw("Same");
drawText("anti-k_{T} PF R=0.3",0.2,0.90);
drawText("|#eta_{barrel}|<1.3",0.2,0.84);
drawText("|#Delta#phi_{1,2}|>2.5",0.2,0.77);
drawText(Form("%d<p_{T,ave}<%d GeV/c",ptlow,pthigh),0.2,0.71);
c3->SaveAs(Form("plots/ratios_%s_%s_%s_%s.png",mode,binning,algo,dataset));
c3->SaveAs(Form("plots/ratios_%s_%s_%s_%s.pdf",mode,binning,algo,dataset));
TFile *outf = new TFile(Form("Corrections/Casym_%s_%s_algo_%s_pt%d_%d_%s.root",mode,binning,algo,ptlow,pthigh,dataset),"recreate");
C_asym->Write();
R_data->Write();
R_mc->Write();
outf->Close();
}