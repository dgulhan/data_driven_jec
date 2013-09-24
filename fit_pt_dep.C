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
 
void pt_dependence(){
TH1D::SetDefaultSumw2();
int npt=6;
int ptlow[]={20,50,70,90,130,200}; 
int pthigh[]={50,70,90,130,200,300};
char *mode="pp";
char *algo = "ak3PF";
char *binning ="9bins";
TString dataset[]={"jet20","jet40","jet40","jet80","jet80","jet80"};
TFile *file[npt];
int linesty[]={2,1,1,1,1,1};
// int col[]={1,kBlue,kAzure+6,kGreen+1,kOrange-3,kRed};
int col[]={kBlue,kAzure+6,kGreen+1,kOrange-3,kRed,kOrange-3,kGreen+1,kAzure+6,kBlue};
int sty[] = {24,25,26,27,29,33,22,21,20};
int size[] = {1,1,1,2,2,2,1,1,1};
 
TLegend *t8=new TLegend(0.55,0.6,0.94,0.93); 
t8->SetFillColor(0);
t8->SetBorderSize(0);
t8->SetFillStyle(0);
t8->SetTextFont(63);
t8->SetTextSize(18); 
int netabins=9;
double etabins_9bins[] = {-3, -2,   -1.3,   -0.8,  -0.4,   0.4,     0.8,  1.3,     2,    3};

double alphacut_high=0.2;
TH1D *C_asym[npt];
TH1D *h_pt_ave[npt];
double pt[npt];
double pterr[npt];
double corr[netabins][npt];
double correrr[netabins][npt];
TGraphErrors *corr_vs_pt[netabins];
for(int ipt=0;ipt<npt;ipt++){
 file[ipt]=new TFile(Form("Corrections/Casym_%s_%s_algo_%s_pt%d_%d_%s_alphahigh_%d_phicut250.root",mode,binning,algo,ptlow[ipt],pthigh[ipt],dataset[ipt].Data(),(int)(alphacut_high*100)));
 C_asym[ipt]=(TH1D*)file[ipt]->Get("C_asym");
 h_pt_ave[ipt]=(TH1D*)file[ipt]->Get("h_pt_ave");
 pt[ipt]=h_pt_ave[ipt]->GetMean();
 pterr[ipt]=h_pt_ave[ipt]->GetMeanError();
 for(int ieta=0;ieta<netabins;ieta++){
  corr[ieta][ipt] = C_asym[ipt]->GetBinContent(ieta+1);
  correrr[ieta][ipt] = C_asym[ipt]->GetBinError(ieta+1);
 }
}

for(int ieta=0;ieta<netabins;ieta++){
 corr_vs_pt[ieta]=new TGraphErrors(npt,pt,corr[ieta],pterr,correrr[ieta]);
 corr_vs_pt[ieta]->SetMarkerColor(col[ieta]);
 corr_vs_pt[ieta]->SetLineColor(col[ieta]);
 corr_vs_pt[ieta]->SetMarkerStyle(sty[ieta]);
 corr_vs_pt[ieta]->SetMarkerSize(size[ieta]);
 if(ieta>3)t8->AddEntry(corr_vs_pt[ieta],Form("%.1f<|eta|<%.1f",fabs(etabins_9bins[ieta]),fabs(etabins_9bins[ieta+1])),"p");
}

// TCanvas *c1 = new TCanvas("c1","",600,600);
// TH1D * empty = new TH1D("empty","",10,20,300);
// empty->Fill(100,1000);
// empty->SetMaximum(1.2);
// empty->SetMinimum(0.95);
// empty->GetXaxis()->SetTitle("p_{T}^{ave}");
// empty->GetXaxis()->CenterTitle();
// empty->GetYaxis()->SetTitle("R_{MC}/R_{data}");
// empty->GetYaxis()->CenterTitle();
// empty->Draw();
// for(int ieta=0;ieta<netabins;ieta++){
// corr_vs_pt[ieta]->Draw("same p");
// }
// t8->Draw("same");
 // drawText("pp #sqrt{s}=2.76 TeV",0.23,0.89);
 // drawText("anti-k_{T}, R=0.3",0.23,0.81);
 // drawText("no UE subt",0.23,0.73);
// c1->SaveAs("fit_pt_dep.png");

TF1 * fit[netabins];
for(int ieta=0;ieta<(netabins);ieta++){
 if(ieta==0 || ieta==(netabins-1)){
  fit[ieta]=new TF1(Form("netabins_%d",ieta),"[0]*(x-[1])*(x-[2])",20,300);
  fit[ieta]->SetParameters(-0.001,0,200);
 }
 else{
 fit[ieta]=new TF1(Form("netabins_%d",ieta),"[0]*x+[1]",20,300);
 fit[ieta]->SetParameters(0.1,1);
 }
 corr_vs_pt[ieta]->Fit(fit[ieta],"R LL");
}
for(int ieta=0; ieta<netabins ; ieta++){
 if(ieta==0 || ieta==(netabins-1)) cout<<etabins_9bins[ieta]<<"-"<<etabins_9bins[ieta+1]<< " p0 = "<<fit[ieta]->GetParameter(0)<<" p1 = "<<fit[ieta]->GetParameter(1)<<" p2 = "<<fit[ieta]->GetParameter(2)<<endl;
 else cout<<etabins_9bins[ieta]<<"-"<<etabins_9bins[ieta+1]<<" slope = "<<fit[ieta]->GetParameter(0)<<" y-intercept = "<<fit[ieta]->GetParameter(1)<<endl;
}

TCanvas *c2 = new TCanvas("c2","",600,600);
corr_vs_pt[0]->Draw("Ap");
c2->SaveAs("fitcheck0.png");
TCanvas *c3 = new TCanvas("c3","",600,600);
corr_vs_pt[netabins-1]->Draw("Ap");
c3->SaveAs("fitcheck8.png");
}