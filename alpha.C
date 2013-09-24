#include "relativeResponse_MC.C"

void alpha(){
TH1D::SetDefaultSumw2();
double etacut=3;
double phicut=2.5;
int docorrection[] ={1,0};
char *mode="pPb";//"pPb","pp","Pbp";
char *algo="akPu3PF";
TString infile;
double alphacut_low=0;
double alphacut_high=0.2;
int ndataset=2;
TString dataset[]={"MB","jet20","jet80"};
// TString dataset[]={"jet40","jet80"};
int ptlow = 20;
int pthigh = 40;

// double etalow=0.0;
// double etahigh=1.3;
double etalow=1.3;
double etahigh=3;

TFile *f_data[ndataset];
TFile *f_mc = new TFile(Form("ntuples/ntuple_relativeResponse_MC_%s_eta%d_algo%s.root",mode,(int)etacut,algo));
for(int idataset=0; idataset<ndataset; idataset++){
 f_data[idataset] = new TFile(Form("ntuples/ntuple_relativeResponse_%s_eta%d_corrected%d_%s_%s.root",mode,(int)etacut,docorrection[idataset],dataset[idataset].Data(),algo));
}

int col[3] = {kBlue,kGreen+1,kRed};
// int col[3] = {kGreen+1,kRed};
int sty[5] = {24,25,26,28,30};

 
 TLegend *t8=new TLegend(0.55,0.7,0.94,0.93); 
 t8->SetFillColor(0);
 t8->SetBorderSize(0);
 t8->SetFillStyle(0);
 t8->SetTextFont(63);
 t8->SetTextSize(18); 
  
 TLegend *t9=new TLegend(0.55,0.7,0.94,0.93); 
 t9->SetFillColor(0);
 t9->SetBorderSize(0);
 t9->SetFillStyle(0);
 t9->SetTextFont(63);
 t9->SetTextSize(18); 

TTree *nt_mc = (TTree*)f_mc->Get("ntdijet");
TTree *nt_data[ndataset];
TH1D * alpha_data[ndataset];
TH1D * B_data[ndataset];
TH1D * B_mc = new TH1D("B_mc","",11,-1,1);
nt_mc->Draw("(ptProbe-ptBarrel)/ptAverage>>B_mc",Form("weight*(dphi>2.5 && ptAverage>=%d && ptAverage<%d && abs(etaProbe)>%.2f)",ptlow,pthigh,etalow));
B_mc->Scale(1/B_mc->Integral());
for(int idataset=0;idataset<ndataset;idataset++){ 
 nt_data[idataset]= (TTree*)f_data[idataset]->Get("ntdijet");
 alpha_data[idataset] = new TH1D(Form("alpha_%s",dataset[idataset].Data()),"",10,0,1);
 B_data[idataset] = new TH1D(Form("B_%s",dataset[idataset].Data()),"",11,-1,1);
 nt_data[idataset]->Draw(Form("alpha>>alpha_%s",dataset[idataset].Data()),Form("dphi>2.5 && jtpt3>0 && ptAverage>=%d && ptAverage<%d && abs(etaProbe)>%.2f",ptlow,pthigh,etalow));
 // nt_data[idataset]->Draw(Form("alpha>>alpha_%s",dataset[idataset].Data()),Form("dphi>2.5 && jtpt3>0 && ptAverage>=%d && ptAverage<%d && abs(etaProbe)<%.2f",ptlow,pthigh,etahigh));
 nt_data[idataset]->Draw(Form("(ptProbe-ptBarrel)/ptAverage>>B_%s",dataset[idataset].Data()),Form("dphi>2.5 && ptAverage>=%d && ptAverage<%d && abs(etaProbe)>%.2f",ptlow,pthigh,etalow));
 // nt_data[idataset]->Draw(Form("(ptProbe-ptBarrel)/ptAverage>>B_%s",dataset[idataset].Data()),Form("dphi>2.5 && ptAverage>=%d && ptAverage<%d && abs(etaProbe)<%.2f",ptlow,pthigh,etahigh));
 B_data[idataset]->Scale(1/B_data[idataset]->Integral());
 // B_data[idataset]->SetMarkerSize(0);
 // alpha_data[idataset]->SetMarkerSize(0);
 B_data[idataset]->SetMarkerColor(col[idataset]);
 B_data[idataset]->SetLineColor(col[idataset]);
 B_data[idataset]->SetMarkerStyle(sty[idataset]);
 alpha_data[idataset]->SetMarkerColor(col[idataset]);
 alpha_data[idataset]->SetMarkerStyle(sty[idataset]);
 alpha_data[idataset]->SetLineColor(col[idataset]);
 alpha_data[idataset]->Scale(1/alpha_data[idataset]->Integral());
 t8->AddEntry(B_data[idataset],Form("%s",dataset[idataset].Data()),"p");
 t9->AddEntry(B_data[idataset],Form("%s",dataset[idataset].Data()),"p");
}
B_mc->SetMarkerSize(0);
B_mc->SetFillColor(17);
if(mode=="pPb" || mode=="Pbp") t8->AddEntry(B_mc,"PYTHIA+HIJING","fl");
if(mode=="pp") t8->AddEntry(B_mc,"PYTHIA","fl");

TCanvas * c1 = new TCanvas("c1","",600,600);
B_mc->GetXaxis()->SetTitle("Dijet Balance");
B_mc->GetXaxis()->CenterTitle();
B_mc->GetYaxis()->SetTitle("EventFraction");
B_mc->GetYaxis()->CenterTitle();
B_mc->SetMaximum(0.45);
B_mc->SetMinimum(0.);
B_mc->Draw();
B_mc->Draw("same hist");
for(int idataset=0; idataset<ndataset;idataset++){
 B_data[idataset]->Draw("same");
 // B_data[idataset]->Draw("hist same");
}
t8->Draw("same");
if(mode=="pp") drawText(Form("%s, #sqrt{s_{NN}} = 2.76 TeV",mode),0.25,0.9);
if(mode=="pPb") drawText(Form("%s, #sqrt{s_{NN}} = 5.02 TeV",mode),0.25,0.9);
drawText("anti-k_{T} PF, R=0.3",0.25,0.85);
if(algo=="akPu3PF") drawText("UE subt",0.25,0.8);
if(algo=="ak3PF") drawText("no UE subt",0.25,0.8);
drawText("pp tracking",0.25,0.75);
 drawText(Form("%d #leq p_{T}^{ave} < %d",ptlow,pthigh),0.25,0.7);
 drawText(Form("|#eta_{probe}|>%.2f",etalow),0.25,0.65);
 // drawText(Form("|#eta_{probe}|<%.2f",etahigh),0.25,0.65);

c1->RedrawAxis();
c1->SaveAs(Form("B_trigger_bias_%s_%s_pt%d_%d_eta_%d_%d.png",mode,algo,ptlow,pthigh,(int)(10*etalow),(int)(etahigh*10)));
c1->SaveAs(Form("B_trigger_bias_%s_%s_pt%d_%d_eta_%d_%d.pdf",mode,algo,ptlow,pthigh,(int)(10*etalow),(int)(etahigh*10)));

TCanvas *c2 = new TCanvas("c2","",600,600);
alpha_data[0]->SetMaximum(0.45);
alpha_data[0]->SetMinimum(0);
alpha_data[0]->GetXaxis()->SetTitle("#alpha=p_{T,3}/p_{T}^{ave}");
alpha_data[0]->GetXaxis()->CenterTitle();
alpha_data[0]->GetYaxis()->SetTitle("Event Fraction");
alpha_data[0]->GetYaxis()->CenterTitle();
alpha_data[0]->Draw();
// alpha_data[0]->Draw("same hist");

for(int idataset=1;idataset<ndataset;idataset++){
 alpha_data[idataset]->Draw("same");
 // alpha_data[idataset]->Draw("same hist");
}
if(mode=="pp") drawText(Form("%s, #sqrt{s_{NN}} = 2.76 TeV",mode),0.25,0.9);
if(mode=="pPb") drawText(Form("%s, #sqrt{s_{NN}} = 5.02 TeV",mode),0.25,0.9);
drawText("anti-k_{T} PF, R=0.3",0.25,0.85);
if(algo=="akPu3PF") drawText("UE subt",0.25,0.8);
if(algo=="ak3PF") drawText("no UE subt",0.25,0.8);
drawText("pp tracking",0.25,0.75);
 drawText(Form("%d #leq p_{T}^{ave} < %d",ptlow,pthigh),0.25,0.7);
 // drawText(Form("|#eta_{probe}|>%.2f",etalow),0.25,0.65);
 drawText(Form("|#eta_{probe}|<%.2f",etahigh),0.25,0.65);

t9->Draw("same");
c2->RedrawAxis();
c2->SaveAs(Form("alpha_trigger_bias_%s_%s_pt%d_%d_eta_%d_%d.png",mode,algo,ptlow,pthigh,(int)(10*etalow),(int)(etahigh*10)));
c2->SaveAs(Form("alpha_trigger_bias_%s_%s_pt%d_%d_eta_%d_%d.pdf",mode,algo,ptlow,pthigh,(int)(10*etalow),(int)(etahigh*10)));

}