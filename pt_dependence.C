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
 
void pt_dependence(){
int npt=6;
int ptlow[]={20,40,60,80,100,140}; 
int pthigh[]={40,60,80,100,140,300};
char *mode="pPb";
char *algo = "akPu3PF";
char *binning ="double_hcalbins";
TString dataset[]={"jet20","jet40","jet40","jet80","jet80","jet80"};
TFile *file[npt];
int linesty[]={2,1,1,1,1,1};
int col[]={1,kBlue,kAzure+6,kGreen+1,kOrange-3,kRed};
int sty[] = {22,22,22,24,25,26};
 
 TLegend *t8=new TLegend(0.55,0.6,0.94,0.93); 
 t8->SetFillColor(0);
 t8->SetBorderSize(0);
 t8->SetFillStyle(0);
 t8->SetTextFont(63);
 t8->SetTextSize(18); 
 
double alphacut_high=0.2;
TH1D *C_asym[npt];
for(int ipt=0;ipt<npt;ipt++){
 file[ipt]=new TFile(Form("Corrections/Casym_%s_%s_algo_%s_pt%d_%d_%s_alphahigh_%d.root",mode,binning,algo,ptlow[ipt],pthigh[ipt],dataset[ipt].Data(),(int)(alphacut_high*100)));
 C_asym[ipt]=(TH1D*)file[ipt]->Get("C_asym");
 C_asym[ipt]->SetLineStyle(linesty[ipt]);
 C_asym[ipt]->SetLineColor(col[ipt]);
 C_asym[ipt]->SetMarkerColor(col[ipt]);
 C_asym[ipt]->SetMarkerStyle(sty[ipt]);
 if(ipt==0) t8->AddEntry(C_asym[ipt],Form("%d #leq p_{T}^{ave} < %d, %s",ptlow[ipt],pthigh[ipt],dataset[ipt].Data()),"l");
 if(ipt==1 || ipt==2) t8->AddEntry(C_asym[ipt],Form("%d #leq p_{T}^{ave} < %d, %s",ptlow[ipt],pthigh[ipt],dataset[ipt].Data()),"l");
 if(ipt>2) t8->AddEntry(C_asym[ipt],Form("%d #leq p_{T}^{ave} < %d, %s",ptlow[ipt],pthigh[ipt],dataset[ipt].Data()),"p");
}
C_asym[0]->SetMarkerSize(0);
C_asym[1]->SetMarkerSize(0);
C_asym[2]->SetMarkerSize(0);

TCanvas *c1 = new TCanvas("c1","",600,600);
C_asym[0]->SetMaximum(1.4);
C_asym[0]->SetMinimum(0.9);
C_asym[0]->Draw();
C_asym[0]->Draw("same hist");
for(int ipt=1;ipt<npt;ipt++){
 if(ipt<3){
  C_asym[ipt]->Draw("same");
  C_asym[ipt]->Draw("same hist");
 }else{
  C_asym[ipt]->Draw("same");
 }
}
 if(mode=="pp") drawText(Form("%s, #sqrt{s_{NN}} = 2.76 TeV",mode),0.25,0.9);
 if(mode=="pPb") drawText(Form("%s, #sqrt{s_{NN}} = 5.02 TeV",mode),0.25,0.9);
 drawText("anti-k_{T} PF, R=0.3",0.25,0.85);
 if(algo=="akPu3PF") drawText("UE subt",0.25,0.8);
 if(algo=="ak3PF") drawText("no UE subt",0.25,0.8);
 drawText("pp tracking",0.25,0.75);
t8->Draw("same");
c1->SaveAs(Form("pt_dependence_%s.png",mode));
c1->SaveAs(Form("pt_dependence_%s.pdf",mode));
}