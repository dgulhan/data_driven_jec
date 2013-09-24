#include <iostream>
#include <vector>
#include <algorithm>

#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom1.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include "TNtuple.h"
#include "TLine.h"
#include "TF1.h"
 
#include "/d102/dgulhan/PA2013Analysis/Dijet/d0127/plotFigure.C"
 
void calculateCasym_vs_pt(){
TH1D::SetDefaultSumw2();
int npt=6;
int ptlow[]={20,50,70,90,130,200}; 
int pthigh[]={50,70,90,130,200,300};
char *mode="pPb";
char *algo = "ak3PF";
int docorrection=1;
char * dataset[]={"jet20","jet40","jet40","jet80","jet80","jet80"};
TFile *f_data[npt];
double pt[npt];
double pterr[npt];
double corr[npt];
double correrr[npt];
double etacut=3;
TFile * f_mc = new TFile(Form("ntuples/ntuple_relativeResponse_MC_%s_eta%d_algo%s.root",mode,(int)etacut,algo));
for(int ipt=0;ipt<npt;ipt++){
 f_data[ipt] = new TFile(Form("ntuples/ntuple_relativeResponse_%s_eta%d_corrected%d_%s_%s.root",mode,(int)etacut,docorrection,dataset[ipt],algo));
}
TTree * nt_mc=(TTree*)f_mc->Get("ntdijet");

TTree *nt_data[npt];
TH1D * B_data[npt];
TH1D * B_mc[npt];
TH1D * h_pt[npt];
double alphacut_low=0;
double alphacut_high=0.2;
double R_data[npt];
double Rerr_data[npt];
double R_mc[npt];
double Rerr_mc[npt];
double phicut=2.5;


for(int ipt=0;ipt<npt;ipt++){
 // nt_data[ipt]=(TTree*)f_data[ipt]->Get("ntdijet_corr");
 nt_data[ipt]=(TTree*)f_data[ipt]->Get("ntdijet_corr");
 B_data[ipt]=new TH1D(Form("B_data_%d",ipt),"",60,-3,3);
 B_mc[ipt]=new TH1D(Form("B_mc_%d",ipt),"",60,-3,3);
 h_pt[ipt]=new TH1D(Form("h_pt_%d",ipt),"",40,ptlow[ipt],pthigh[ipt]);
 
 nt_mc->Draw(Form("(ptProbe-ptBarrel)/ptAverage>>B_mc_%d",ipt),Form("weight*(abs(etaProbe)<2 && ptAverage>%d && ptAverage<%d && ( (jtpt3>0 && (alpha>=%.2f && alpha<%.2f)) || (jtpt3<0)) && dphi>2.5 && pthat>%d && (ptAverage<80 || weight<0.0708746))",ptlow[ipt],pthigh[ipt],alphacut_low,alphacut_high,(ptlow-10)));
 // nt_mc->Draw(Form("(ptProbe-ptBarrel)/ptAverage>>B_mc_%d",ipt),Form("weight*(abs(etaProbe)<2 && ptAverage>%d && ptAverage<%d && ( (jtpt3>0 && (alpha>=%.2f && alpha<%.2f)) || (jtpt3<0)) && dphi>2.5 && pthat>%d && weight<0.0037981)",ptlow[ipt],pthigh[ipt],alphacut_low,alphacut_high,(ptlow-10)));
 nt_data[ipt]->Draw(Form("(ptProbe-ptBarrel)/ptAverage>>B_data_%d",ipt),Form("abs(etaProbe)<2 && ptAverage>%d && ptAverage<%d && ((jtpt3>0 && (alpha>=%.2f && alpha<%.2f)) || (jtpt3<0)) && dphi>%.2f",ptlow[ipt],pthigh[ipt],alphacut_low,alphacut_high,phicut));
 nt_data[ipt]->Draw(Form("(ptProbe+ptBarrel)/2>>h_pt_%d",ipt),Form("abs(etaProbe)<2 && ptAverage>%d && ptAverage<%d && ((jtpt3>0 && (alpha>=%.2f && alpha<%.2f)) || (jtpt3<0)) && dphi>%.2f",ptlow[ipt],pthigh[ipt],alphacut_low,alphacut_high,phicut));
 pt[ipt]=h_pt[ipt]->GetMean();
 pterr[ipt]=h_pt[ipt]->GetMeanError();
 B_mc[ipt]->Scale(1/B_mc[ipt]->Integral());
 R_data[ipt]=(2+B_data[ipt]->GetMean())/(2-B_data[ipt]->GetMean());
 R_mc[ipt]=(2+B_mc[ipt]->GetMean())/(2-B_mc[ipt]->GetMean());
 
 Rerr_mc[ipt] = 4*B_mc[ipt]->GetMeanError()/pow((2-B_mc[ipt]->GetMean()),2);
 Rerr_data[ipt] = 4*B_data[ipt]->GetMeanError()/pow((2-B_data[ipt]->GetMean()),2);
 cout << ipt<<"mc error "<<Rerr_mc[ipt]<<" data error "<<Rerr_data[ipt]<<endl;
 cout << ipt<<"mc B error "<<B_mc[ipt]->GetMeanError()<<" mc B mean "<<(double)B_mc[ipt]->GetMean()<<endl;
 cout << ipt<<"data B error "<<B_data[ipt]->GetMeanError()<<" data B mean "<<(double)B_data[ipt]->GetMean()<<endl;
 corr[ipt]=R_mc[ipt]/R_data[ipt];
 correrr[ipt]=corr[ipt]*sqrt(pow(Rerr_mc[ipt]/R_mc[ipt],2)+pow(Rerr_data[ipt]/R_data[ipt],2));
 
}
TGraphErrors * g = new TGraphErrors(npt,pt,corr,pterr,correrr);

TCanvas *c1 = new TCanvas("c1","",600,600);
TH1D * empty = new TH1D("empty","",10,20,300);
empty->Fill(100,1000);
empty->SetMaximum(1.05);
empty->SetMinimum(0.95);
empty->Draw();
empty->GetXaxis()->SetTitle("p_{T}^{ave}");
empty->GetXaxis()->CenterTitle();
empty->GetYaxis()->SetTitle("R_{MC}/R^{corr}_{data}");
empty->GetYaxis()->CenterTitle();
g->Draw("same p");
c1->SaveAs(Form("plots_JEC/Step2_Corr_%s_%s.png",mode,algo));

TF1* f = new TF1("f","1-[0]/pow(x,[1])",20,200);
f->SetParameters(0.9,0.01);
g->Fit(f,"R");

TCanvas *c2 = new TCanvas("c2","",600,600);
empty->Draw();
g->Draw("same p");
c2->SaveAs(Form("plots_JEC/Step2_Corr_fit_%s_%s.png",mode,algo));

TCanvas *c3 = new TCanvas("c3","",600,600);
B_mc[4]->SetMarkerColor(2);
B_mc[4]->SetLineColor(2);
B_mc[5]->Draw();

B_mc[4]->Draw("same");
c3->SaveAs("plots_JEC/B4.png");
}