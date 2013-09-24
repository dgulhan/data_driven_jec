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
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TLine.h"
#include "TF1.h"
 
 #include "/d102/dgulhan/PA2013Analysis/Dijet/d0127/plotFigure.C"
 
void pt_dependence(){
TH1D::SetDefaultSumw2();
int npt=12;
// int ptlow[]={20,30,40,50,60,70,80,100,120,140,180,250};
// int pthigh[]={30,40,50,60,70,80,100,120,140,180,250,300};
int ptlow[]={25,45,50,60,70,80,100,120,140,180,250};
int pthigh[]={45,50,60,70,80,100,120,140,180,250,300};
char *mode ="pPb";
char *algo;
int doalphacut=1;
int docorrection=1;
char *MC;
if(mode=="pPb"){
 algo ="ak3PF";
 MC="PYTHIA+HIJING";
}
if(mode=="pp"){ 
 algo ="ak3PF";
 MC="PYTHIA";
}
char *dataset[3]={"jet20","jet40","jet80"};
TFile *file[npt];
int linesty[]={2,1,1,1,1,1};
int col[]={1,kBlue,kAzure+6,kGreen+1,kOrange-3,kRed};
int sty[] = {22,22,22,24,25,26};

double sigma_data[npt];
double sigma_mc[npt];
double sigma_data_err[npt];
double sigma_mc_err[npt];
double pt_ave[npt];
double pt_ave_err[npt];
TLegend *t8=new TLegend(0.55,0.7,0.94,0.93); 
t8->SetFillColor(0);
t8->SetBorderSize(0);
t8->SetFillStyle(0);
t8->SetTextFont(63);
t8->SetTextSize(18); 
 double diff[npt];
 double diff_err[npt];
 double diff_quad[npt];
 double diff_quad_err[npt];
double alphacut_high=0.2;
double alphacut_low=0;
TFile *f_mc = new TFile(Form("ntuples/ntuple_relativeResponse_MC_%s_eta3_algo%s.root",mode,algo));
TFile *f_data[3];
TTree *nt_data[3];
for(int i=0;i<3;i++){
 f_data[i] = new TFile(Form("ntuples/ntuple_relativeResponse_%s_eta3_corrected%d_%s_%s.root",mode,docorrection,dataset[i],algo));
 nt_data[i] = (TTree*)f_data[i]->Get("ntdijet");
}
TTree *nt_mc = (TTree*)f_mc->Get("ntdijet");
TH1D *B_data[npt];
TH1D *B_mc[npt];
TH1D *hpt_data[npt];
TH1D *hpt_mc[npt];
 

for(int ipt=0;ipt<npt;ipt++){
 B_data[ipt] = new TH1D(Form("B_data_%d",ipt),"",100,-3,3);
 B_mc[ipt] = new TH1D(Form("B_mc_%d",ipt),"",100,-3,3);
 hpt_data[ipt] = new TH1D(Form("hpt_data_%d",ipt),"",50,ptlow[ipt],pthigh[ipt]);
 hpt_mc[ipt] = new TH1D(Form("hpt_mc_%d",ipt),"",50,ptlow[ipt],pthigh[ipt]);
 int i;

 if(ipt<1) i=0;
 else if(ipt<5) i=1;
 else i=2;
 
 if(doalphacut) nt_data[i]->Draw(Form("(ptProbe-ptBarrel)/ptAverage>>B_data_%d",ipt),Form("dphi>2.5 && ptAverage>=%d && ptAverage<%d && ((jtpt3>0 && (alpha>=%.2f && alpha<%.2f)) || (jtpt3<0)) && abs(etaProbe)<1.3",ptlow[ipt],pthigh[ipt],alphacut_low,alphacut_high));
 if(!doalphacut) nt_data[i]->Draw(Form("(ptProbe-ptBarrel)/ptAverage>>B_data_%d",ipt),Form("dphi>2.5 && ptAverage>=%d && ptAverage<%d && abs(etaProbe)<1.3",ptlow[ipt],pthigh[ipt]));
 if(doalphacut) nt_data[i]->Draw(Form("ptAverage>>hpt_data_%d",ipt),Form("dphi>2.5 && ptAverage>=%d && ptAverage<%d && ((jtpt3>0 && (alpha>=%.2f && alpha<%.2f)) || (jtpt3<0))&& abs(etaProbe)<1.3",ptlow[ipt],pthigh[ipt],alphacut_low,alphacut_high));
 if(!doalphacut) nt_data[i]->Draw(Form("ptAverage>>hpt_data_%d",ipt),Form("dphi>2.5 && ptAverage>=%d && ptAverage<%d && abs(etaProbe)<1.3",ptlow[ipt],pthigh[ipt]));
 if(doalphacut)nt_mc->Draw(Form("(ptProbe-ptBarrel)/ptAverage>>B_mc_%d",ipt),Form("weight*(dphi>2.5 && ptAverage>=%d && ptAverage<%d && ((jtpt3>0 && (alpha>=%.2f && alpha<%.2f)) || (jtpt3<0)) && pthat>%d && abs(etaProbe)<1.3)",ptlow[ipt],pthigh[ipt],alphacut_low,alphacut_high,(ptlow[ipt]-20)));
 if(doalphacut)nt_mc->Draw(Form("ptAverage>>hpt_mc_%d",ipt),Form("weight*(dphi>2.5 && ptAverage>=%d && ptAverage<%d && ((jtpt3>0 && (alpha>=%.2f && alpha<%.2f)) || (jtpt3<0)) && pthat>%d && abs(etaProbe)<1.3)",ptlow[ipt],pthigh[ipt],alphacut_low,alphacut_high,(ptlow[ipt]-20)));
 if(!doalphacut) nt_mc->Draw(Form("(ptProbe-ptBarrel)/ptAverage>>B_mc_%d",ipt),Form("weight*(dphi>2.5 && ptAverage>=%d && ptAverage<%d && pthat>%d && abs(etaProbe)<1.3)",ptlow[ipt],pthigh[ipt],(ptlow[ipt]-20)));
 if(!doalphacut)nt_mc->Draw(Form("ptAverage>>hpt_mc_%d",ipt),Form("weight*(dphi>2.5 && ptAverage>=%d && ptAverage<%d && pthat>%d && abs(etaProbe)<1.3)",ptlow[ipt],pthigh[ipt],(ptlow[ipt]-20)));
 B_data[ipt]->Scale(1/B_data[ipt]->Integral());
 B_mc[ipt]->Scale(1/B_mc[ipt]->Integral());
 sigma_data[ipt]=B_data[ipt]->GetRMS(); 
 sigma_mc[ipt]=B_mc[ipt]->GetRMS(); 
 sigma_data_err[ipt]=B_data[ipt]->GetRMSError(); 
 sigma_mc_err[ipt]=B_mc[ipt]->GetRMSError(); 
 pt_ave[ipt]=hpt_data[ipt]->GetMean();
 pt_ave_err[ipt]=hpt_data[ipt]->GetMeanError();
 diff[ipt]=sigma_data[ipt]-sigma_mc[ipt];
 diff_quad[ipt]=sqrt(pow(sigma_data[ipt],2)-pow(sigma_mc[ipt],2))/sqrt(2);
 diff_err[ipt]=sqrt(pow(sigma_data_err[ipt],2)+pow(sigma_mc_err[ipt],2.));
 diff_quad_err[ipt]=sqrt(pow(sigma_data_err[ipt]*sigma_data[ipt]/diff_quad[ipt],2)+pow(sigma_mc_err[ipt]*sigma_mc[ipt]/diff_quad[ipt],2.))/sqrt(2);
}


TGraphErrors *hsigB_pt_data=new TGraphErrors(npt,pt_ave,sigma_data,pt_ave_err,sigma_data_err);
TGraphErrors *hsigB_pt_mc=new TGraphErrors(npt,pt_ave,sigma_mc,pt_ave_err,sigma_mc_err);
TGraphErrors *hsigB_pt_diff=new TGraphErrors(npt,pt_ave,diff,pt_ave_err,diff_err);
TGraphErrors *hsigB_pt_diff_quad=new TGraphErrors(npt,pt_ave,diff_quad,pt_ave_err,diff_quad_err);

hsigB_pt_data->SetMarkerColor(2);
hsigB_pt_mc->SetMarkerStyle(24);
t8->AddEntry(hsigB_pt_data,mode,"p");
t8->AddEntry(hsigB_pt_mc,MC,"p");

TH1D * empty = new TH1D("empty","",10,20,300);
empty->Fill(100,1000);
empty->SetMaximum(0.5);
empty->SetMinimum(0);

TCanvas *c1 = new TCanvas("c1","",600,600);
empty->GetXaxis()->SetTitle("p_{T}^{ave}");
empty->GetYaxis()->SetTitle("#sigma(B)");
empty->Draw();
hsigB_pt_data->Draw("same p");
hsigB_pt_mc->Draw("same p");
t8->Draw("same");
c1->SaveAs(Form("plots_JER/resol_pt_dependence_%s_alphacut%d_%s.png",mode,doalphacut,algo));
// c1->SaveAs(Form("plots_JER/resol_pt_dependence_%s_alphacut%d_%s_onlyjet20.png",mode,doalphacut,algo));

TCanvas *c2 = new TCanvas("c2","",600,600);
empty->SetMaximum(0.5);
empty->GetXaxis()->SetTitle("p_{T}^{ave}");
empty->GetYaxis()->SetTitle("#sigma(B_{data})-#sigma(B_{mc})");
empty->Draw();
hsigB_pt_diff->Draw("same p");
c2->SaveAs(Form("plots_JER/resol_pt_dependence_diff_%s_alphacut%d_%s.png",mode,doalphacut,algo));
// c2->SaveAs(Form("plots_JER/resol_pt_dependence_diff_%s_alphacut%d_%s_onlyjet20.png",mode,doalphacut,algo));

TCanvas *c2_2 = new TCanvas("c2_2","",600,600);
empty->SetMaximum(0.3);
empty->GetXaxis()->SetTitle("p_{T}^{ave}");
empty->GetYaxis()->SetTitle("#sqrt{#sigma(B_{data})^{2}-#sigma(B_{mc})^{2}}/#sqrt{2}");
empty->Draw();
hsigB_pt_diff_quad->Draw("same p");
// drawText(Form("%s trigger",dataset),0.22,0.9);
drawText(Form("#Delta#phi > 2.5"),0.22,0.83);
drawText(Form("#alpha < 0.2"),0.22,0.76);
c2_2->SaveAs(Form("plots_JER/resol_pt_dependence_diff_quad_%s_alphacut%d_%s.png",mode,doalphacut,algo));
// c2_2->SaveAs(Form("plots_JER/resol_pt_dependence_diff_quad_%s_alphacut%d_%s_onlyjet20.png",mode,doalphacut,algo));

TCanvas *c3[npt];
for(int ipt=0;ipt<npt;ipt++){
 int i;
 if(ipt<1) i=0;
 else if(ipt<5) i=1;
 else i=2;
 
 c3[ipt]=new TCanvas(Form("c3_%d",ipt),"",600,600);
 B_data[ipt]->Draw();
 B_mc[ipt]->Draw("same hist");
 c3[ipt]->SaveAs(Form("plots_JER/B_%d_%s_alphacut%d_%s.png",ipt,mode,doalphacut,dataset[i]));
}

// TH2D *h_phi_data = new TH2D("h_phi_data","",40,2.5,TMath::Pi(),40,0,1);
// TH2D *h_phi_mc = new TH2D("h_phi_mc","",40,2.5,TMath::Pi(),40,0,1);
// nt_data[i]->Draw("alpha:dphi>>h_phi_data","jtpt3>0 && ptAverage>50");
// nt_mc->Draw("alpha:dphi>>h_phi_mc","weight*(jtpt3>0 && ptAverage>50)");

TFile *outf = new TFile(Form("resol_diff_ptlow25_alphacut%d_%s_%s.root",doalphacut,mode,algo),"recreate");
// TFile *outf = new TFile(Form("resol_diff_ptlow25_alphacut%d_%s_%s_onlyjet20.root",doalphacut,mode,algo),"recreate");
hsigB_pt_data->SetName("data_resol");
hsigB_pt_mc->SetName("mc_resol");
hsigB_pt_diff_quad->SetName("add_smearing");
hsigB_pt_diff_quad->Write();
hsigB_pt_data->Write();
hsigB_pt_mc->Write();
outf->Close();

// TCanvas *c4 = new TCanvas("c4","",600,600);
// h_phi_data->GetXaxis()->SetTitle("#Delta#phi");
// h_phi_data->GetYaxis()->SetTitle("#alpha");
// h_phi_data->Draw("colz");
// c4->SaveAs("h_phi_alpha_data.png");

// TCanvas *c5 = new TCanvas("c5","",600,600);
// h_phi_mc->GetXaxis()->SetTitle("#Delta#phi");
// h_phi_mc->GetYaxis()->SetTitle("#alpha");
// h_phi_mc->Draw("colz");
// c5->SaveAs("h_phi_alpha_mc.png");

TF1 * f = new TF1("f","[0]/pow(x,[1])",50,300);
f->SetParameters(0.9,0.1);
hsigB_pt_diff_quad->Fit(f,"R LL");

TCanvas *c6 = new TCanvas("c6","",600,600);
empty->Draw();
hsigB_pt_diff_quad->Draw("same p");
c6->SaveAs(Form("plots_JER/extra_smearing_fit_%s_%s.png",mode,algo));
// c6->SaveAs(Form("plots_JER/extra_smearing_fit_%s_%s_onlyjet20.png",mode,algo));

}