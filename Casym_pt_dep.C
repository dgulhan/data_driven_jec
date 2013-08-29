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

#include "/d102/dgulhan/PA2013Analysis/Dijet/d0127/plotFigure.C"

void Casym_pt_dep(){
// const int nfile = 4;
// int ptlow[nfile]={40,60,80,100};
// int pthigh[nfile]={60,80,100,140};
const int nfile = 2;
// int ptlow[nfile]={20,30,20,40,40};
// int pthigh[nfile]={30,40,40,300,300};
int ptlow[nfile]={40,40};
int pthigh[nfile]={300,300};
// int ptlow[nfile]={20,40,70,100};
// int pthigh[nfile]={40,70,100,140};

TFile *file[nfile];
TH1D *Casym[nfile];
TH1D *Rdata[nfile];
TH1D *Rmc[nfile];

int col[5] = {kBlue,kAzure+6,kGreen+1,kOrange-3,kRed};
int sty[5] = {24,25,26,28,30};
// TLegend *t8=new TLegend(0.3,0.8,0.94,0.93); 
TLegend *t8=new TLegend(0.45,0.7,0.94,0.93); 
t8->SetFillColor(0);
t8->SetBorderSize(0);
t8->SetFillStyle(0);
t8->SetTextFont(63);
t8->SetTextSize(18); 

for(int ifile=0;ifile<nfile;ifile++){
 // file[ifile] = new TFile(Form("Corrections/Casym_pPb_Aug20_double_hcalbins_pt%d_%d.root",ptlow[ifile],pthigh[ifile]));
 if(ifile<nfile-1)file[ifile] = new TFile(Form("Corrections/Casym_pPb_Aug20_double_hcalbins_pt%d_%d.root",ptlow[ifile],pthigh[ifile]));
 else if(ifile==nfile-1) file[ifile] = new TFile("Corrections/Casym_pPb_Aug20_double_hcalbins_pt40_300_jet40.root");
 //***file[ifile] = new TFile(Form("Corrections/Casym_pp_Aug26_double_hcalbins_pt%d_%d.root",ptlow[ifile],pthigh[ifile]));
 // else if(ifile==2) file[ifile] = new TFile("Corrections/Casym_pPb_Aug20_double_hcalbins_pt70_80_jet60.root");
 // else if(ifile==nfile-1) file[ifile] = new TFile("Corrections/Casym_pPb_Jul17_double_hcalbins.root");
 Casym[ifile] = (TH1D*)file[ifile]->Get("C_asym");
 // Rdata[ifile] = (TH1D*)file[ifile]->Get("R_data");
 // Rmc[ifile] = (TH1D*)file[ifile]->Get("R_mc");
 Casym[ifile]->SetMarkerStyle(sty[ifile]);
 Casym[ifile]->SetMarkerColor(col[ifile]);
 Casym[ifile]->SetLineColor(col[ifile]);
 // Rdata[ifile]->SetMarkerStyle(sty[ifile]);
 // Rdata[ifile]->SetMarkerColor(col[ifile]);
 // Rdata[ifile]->SetLineColor(col[ifile]);
 // Rmc[ifile]->SetMarkerStyle(sty[ifile]);
 // Rmc[ifile]->SetMarkerColor(col[ifile]);
 // Rmc[ifile]->SetLineColor(col[ifile]);
 if(ifile==0)t8->AddEntry(Casym[ifile],Form("%d < p_{T,ave} < %d GeV/c, MB",ptlow[ifile],pthigh[ifile]),"p");
 if(ifile==nfile-1)t8->AddEntry(Casym[ifile],Form("%d < p_{T,ave} < %d GeV/c, jet40",ptlow[ifile],pthigh[ifile]),"p");
 // if(ifile==1 || ifile ==2 )t8->AddEntry(Casym[ifile],Form("%d < p_{T,ave} < %d GeV/c, jet40",ptlow[ifile],pthigh[ifile]),"p");
 // if(ifile==0 || ifile ==3 )t8->AddEntry(Casym[ifile],Form("%d < p_{T,ave} < %d GeV/c, MB",ptlow[ifile],pthigh[ifile]),"p");
 // if(ifile==nfile-1)t8->AddEntry(Casym[ifile],Form("p_{T,ave} > %d GeV/c, jet80",ptlow[ifile]),"p");
}


TCanvas *c1 = new TCanvas("c1","",600,600);
Casym[0]->SetMaximum(1.3);
Casym[0]->Draw();
TLine l(-3,1,3,1);
for(int ifile=1; ifile<nfile; ifile++){
 Casym[ifile]->Draw("same");
}
t8->Draw("same");
l.Draw("same");
drawText("pPb, #sqrt{s_{NN}}=5.02 TeV",0.25,0.9);
drawText("anti-k_{T} PF, R=0.3",0.25,0.85);
drawText("UE subt",0.25,0.8);
drawText("pPb tracking",0.25,0.75);
drawText("Minimum bias data",0.25,0.7);
c1->SaveAs("Casym_pt_dep_pPb_MB.pdf");

// TCanvas *c2 = new TCanvas("c2","",600,600);
// Rdata[0]->Draw();
// for(int ifile=1; ifile<nfile; ifile++){
 // Rdata[ifile]->Draw("same");
// }
// t8->Draw("same");
// c2->SaveAs("Rdata_pt_dep_pPb.pdf");

// TCanvas *c3 = new TCanvas("c3","",600,600);
// Rmc[0]->Draw();
// for(int ifile=1; ifile<nfile; ifile++){
 // Rmc[ifile]->Draw("same");
// }
// t8->Draw("same");
// c3->SaveAs("Rmc_pt_dep_pPb.pdf");
}