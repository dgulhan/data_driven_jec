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
#include "TNtuple.h"
#include "TLine.h"
#include "TF1.h"
 
 #include "/d102/dgulhan/PA2013Analysis/Dijet/d0127/plotFigure.C"

void triggerBias(){
 TH1D::SetDefaultSumw2();
 char *mode="pp";//"pPb","pp","Pbp";
 char *algo="ak3PF";
 // int ndataset=3;
 int ndataset=2;
 int ptlow=40;
 int pthigh=300;
 double alphacut_low=0;
 double alphacut_high=0.2;
 TString dataset[]={"MB","jet40","jet80"};
 // TString dataset[]={"jet40","jet80"};
 // TString dataset[]={"MB","jet40"};
 // char *binning = "double_hcalbins";
 char *binning = "6bins";
 TFile *f[ndataset];
 TH1D* C_asym[ndataset];
 TH1D* R_data[ndataset];
 // int col[5] = {kBlue,kAzure+6,kGreen+1,kOrange-3,kRed};
 int col[3] = {kBlue,kGreen+1,kRed};
 // int col[2] = {kGreen+1,kRed};
 int sty[5] = {24,25,26,28,30};
 
 TLegend *t8=new TLegend(0.45,0.7,0.94,0.93); 
 t8->SetFillColor(0);
 t8->SetBorderSize(0);
 t8->SetFillStyle(0);
 t8->SetTextFont(63);
 t8->SetTextSize(18); 
 
 for(int idataset=0;idataset<ndataset;idataset++){
  f[idataset]=new TFile(Form("Corrections/Casym_%s_%s_algo_%s_pt%d_%d_%s_alphahigh_%d.root",mode,binning,algo,ptlow,pthigh,dataset[idataset].Data(),(int)(alphacut_high*100)));
  C_asym[idataset]=(TH1D*)f[idataset]->Get("C_asym");
  R_data[idataset]=(TH1D*)f[idataset]->Get("R_data");
  for(int i=0;i<R_data[idataset]->GetNbinsX();i++){
   if(R_data[idataset]->GetBinError(i+1)<0.00001){
    C_asym[idataset]->SetBinContent(i+1,1000);
    cout<<1<<endl;
    }
  }
  C_asym[idataset]->SetMarkerColor(col[idataset]);
  C_asym[idataset]->SetLineColor(col[idataset]);
  C_asym[idataset]->SetMarkerStyle(sty[idataset]);
  // C_asym[idataset]->SetMarkerSize(2);
  t8->AddEntry(C_asym[idataset],Form("%s",dataset[idataset].Data()),"p");
 }
 
 TCanvas *c1 = new TCanvas("c1","",600,600);
 C_asym[0]->SetMaximum(1.3);
 C_asym[0]->SetMinimum(0.9);
 C_asym[0]->Draw();
 TLine l(-3,1,3,1);
 for(int ifile=1; ifile<ndataset; ifile++){
  C_asym[ifile]->Draw("same");
 }
 t8->Draw("same");
 l.Draw("same");
 if(mode=="pp") drawText(Form("%s, #sqrt{s_{NN}} = 2.76 TeV",mode),0.25,0.9);
 if(mode=="pPb") drawText(Form("%s, #sqrt{s_{NN}} = 5.02 TeV",mode),0.25,0.9);
 drawText("anti-k_{T} PF, R=0.3",0.25,0.85);
 if(algo=="akPu3PF") drawText("UE subt",0.25,0.8);
 if(algo=="ak3PF") drawText("no UE subt",0.25,0.8);
 drawText("pp tracking",0.25,0.75);
 drawText(Form("%d #leq p_{T}^{ave} < %d",ptlow,pthigh),0.25,0.7);
 c1->SaveAs(Form("plots/triggerBias_%s_%s_pt_%d_%d_alpha_%d_%s.pdf",mode,algo,ptlow,pthigh,(int)(alphacut_high*100),binning));
 c1->SaveAs(Form("plots/triggerBias_%s_%s_pt_%d_%d_alpha_%d_%s.png",mode,algo,ptlow,pthigh,(int)(alphacut_high*100),binning));

}