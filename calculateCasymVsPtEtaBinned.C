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
 
// #include "/d102/dgulhan/PA2013Analysis/Dijet/d0127/plotFigure.C"
 
void calculateCasymVsPtEtaBinned(){
 TH1D::SetDefaultSumw2();
 int npt=7; 
 int ptlow[]= {20,60,80,95,120,150,200}; 
 int pthigh[]={60,80,95,120,150,200,300};
 int etacut=4;
 int netabins=64;
 double etabins[netabins];
 
 double etabins_hcalbins_eta3[]= {-3, -2.853,
                     -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830,
                     -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
                     -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609,
                     -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,  0.000,
                      0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,
                      0.696,  0.783,  0.879,  0.957,  1.044,  1.131,  1.218,
                      1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,
                      1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  2.853,
                      3};
 double etabins_hcalbins_eta4[]= {-4,      -3.664,  -3.314,  -2.964, -2.853,
                     -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830,
                     -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
                     -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609,
                     -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,  0.000,
                      0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,
                      0.696,  0.783,  0.879,  0.957,  1.044,  1.131,  1.218,
                      1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,
                      1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  2.853,
                      2.964,  3.314,  3.664, 4};
 double etabins_double_hcalbins_eta4[]={-4,      -3.664,  -3.314,  -2.964, 
                  -2.650,  -2.322,  -2.043,  -1.830,
                  -1.653,  -1.479, -1.305, 
                  -1.131,  -0.957, -0.783, -0.609,
                   -0.435, -0.261,  -0.087,
                   0.087,  0.261,  0.435,  0.609,
                   0.783,  0.957,  1.131, 
                   1.305,  1.479,  1.653,  1.830,
                   2.043,  2.322,  2.650, 
                   2.964,  3.314,  3.664, 4};
				   
 double etabins_double_hcalbins_eta3[]= {-3,-2.500,  -2.172,-1.740, -1.392, -1.044,  -0.696,  -0.348,  0.000,  0.348,   0.696,   1.044,1.392,   1.740,   2.172,   2.500,   3};
 if(etacut==3){
  netabins=58;
  for(int ieta = 0; ieta<netabins+1; ieta++){
   // etabins_double_hcalbins[ieta] = etabins_double_hcalbins3[ieta];
   etabins[ieta] = etabins_hcalbins_eta3[ieta];
  }
 }
 if(etacut==4){				
  for(int ieta = 0; ieta<netabins+1; ieta++){
   etabins[ieta] = etabins_hcalbins_eta4[ieta];
  }  
 }
 TString mode="pp";
 TString algo = "ak3PF";
 int docorrection=1;
 TString dataset[]={"hcalnoise","jet40 && hcalnoise","jet60 && hcalnoise","jet80 && hcalnoise","jet80 && hcalnoise","jet80 && hcalnoise","jet80 && hcalnoise"};

 int fileindex[] = {0,1,1,2,2,2,2};
 int filemcindex[] = {0,2,2,3,3,4,5};
 int npthat = 6;
 int pthat[] = {15,30,50,80,120,170};

 TFile *f_data[npt];
 double pt[netabins][npt];
 double pterr[netabins][npt];
 double corr[netabins][npt];
 double correrr[netabins][npt];
 TGraphErrors *graph[netabins];
 TGraphErrors *graphclone[netabins];
 TFile * f_mc[npt];  
 
 for(int ipthat = 0; ipthat < npthat; ipthat++){
  f_mc[ipthat] = new TFile(Form("ntuples/pthat%d_%s_eta%d.root",pthat[ipthat],algo.Data(),etacut));
 }

 if(etacut==3){
  f_data[0] = new TFile(Form("ntuples/minbias_%s_eta3_corr.root",algo.Data()));
  f_data[1] = new TFile(Form("ntuples/lowerpt_%s_eta3_corr.root",algo.Data()));
  f_data[2] = new TFile(Form("ntuples/jet80_%s_eta3_corr.root",algo.Data())); 
 }
 if(etacut==4){
  f_data[0] = new TFile(Form("ntuples/minbias_%s_eta4_corr.root",algo.Data()));
  f_data[1] = new TFile(Form("ntuples/lowerpt_%s_eta4_corr.root",algo.Data()));
  f_data[2] = new TFile(Form("ntuples/jet80_%s_eta4_corr.root",algo.Data()));
 }
 TTree * nt_mc[filemcindex[npt-1]+1];
 TTree *nt_data[fileindex[npt-1]+1];

 for(int ifile = 0; ifile < filemcindex[npt-1]+1; ifile++){
  nt_mc[ifile] =(TTree*)f_mc[ifile]->Get("ntdijet");
 }
 for(int ifile = 0; ifile < fileindex[npt-1]+1; ifile++){
  nt_data[ifile] =(TTree*)f_data[ifile]->Get("ntdijet");
 }

 TH1D * B_data[npt][netabins];
 TH1D * B_mc[npt][netabins];
 TH1D * h_pt[npt][netabins];
 double alphacut_low=0;
 double alphacut_high=0.2;
 double R_data[npt][netabins];
 double Rerr_data[npt][netabins];
 double R_mc[npt][netabins];
 double Rerr_mc[npt][netabins];
 double phicut=2.5;
 
 TF1* fit[netabins];
 for(int ieta=0;ieta<netabins;ieta++){
  // fit[ieta] = new TF1(Form("fit%d",ieta)," (1./[2])*([3]*([4]+[5]*[5]*TMath::Log(max([0],min([1],x)))))",20,300);
  // fit[ieta]->SetParameters(55,510,0.997211,1.001024,1.207471,-0.021848);
  fit[ieta] = new TF1(Form("fit%d",ieta)," [0]",20,300);
  fit[ieta]->SetParameter(0,1.);
  for(int ipt=0;ipt<npt;ipt++){
   B_data[ipt][ieta]=new TH1D(Form("B_data_%d_%d",ipt,ieta),"",60,-1.5,1.5);
   B_mc[ipt][ieta]=new TH1D(Form("B_mc_%d_%d",ipt,ieta),"",60,-1.5,1.5);
   h_pt[ipt][ieta]=new TH1D(Form("h_pt_%d_%d",ipt,ieta),"",40,ptlow[ipt],pthigh[ipt]);

   nt_mc[filemcindex[ipt]]->Draw(Form("(ptProbe-ptBarrel)/ptAverage>>B_mc_%d_%d",ipt,ieta),Form("ptAverage>%d && ptAverage<%d && ptProbe>10 && ( (jtpt3>0 && (alpha>=%.2f && alpha<%.2f)) || (jtpt3<0)) && dphi>2.5 && etaProbe>=%.5f && etaProbe<%.5f",ptlow[ipt],pthigh[ipt],alphacut_low,alphacut_high,etabins[ieta],etabins[ieta+1]));

   nt_data[fileindex[ipt]]->Draw(Form("(ptProbe-ptBarrel)/ptAverage>>B_data_%d_%d",ipt,ieta),Form("ptAverage>%d && ptAverage<%d && ptProbe>10 && ((jtpt3>0 && (alpha>=%.2f && alpha<%.2f)) || (jtpt3<0)) && dphi>%.2f&& etaProbe>=%.5f && etaProbe<%.5f&& %s",ptlow[ipt],pthigh[ipt],alphacut_low,alphacut_high,phicut,etabins[ieta],etabins[ieta+1],dataset[ipt].Data()));

   nt_data[fileindex[ipt]]->Draw(Form("(ptProbe+ptBarrel)/2>>h_pt_%d_%d",ipt,ieta),Form("ptAverage>%d && ptAverage<%d && ptProbe>10 && ((jtpt3>0 && (alpha>=%.2f && alpha<%.2f)) || (jtpt3<0)) && dphi>%.2f&& etaProbe>=%.5f && etaProbe<%.5f && %s",ptlow[ipt],pthigh[ipt],alphacut_low,alphacut_high,phicut,etabins[ieta],etabins[ieta+1],dataset[ipt].Data()));

   pt[ieta][ipt]=h_pt[ipt][ieta]->GetMean();
   pterr[ieta][ipt]=h_pt[ipt][ieta]->GetMeanError();
   R_data[ipt][ieta]=(2+B_data[ipt][ieta]->GetMean())/(2-B_data[ipt][ieta]->GetMean());
   R_mc[ipt][ieta]=(2+B_mc[ipt][ieta]->GetMean())/(2-B_mc[ipt][ieta]->GetMean());
   float Berr_mc = sqrt(1./B_mc[ipt][ieta]->Integral());
   float Berr_data = sqrt(1./B_data[ipt][ieta]->Integral());
  
   float Bmean_mc = B_mc[ipt][ieta]->GetMean();
   float Bmean_data = B_data[ipt][ieta]->GetMean();
   Rerr_mc[ipt][ieta] = sqrt(pow(2/(2-Bmean_mc),2.)+pow(2*Bmean_mc/pow(2-Bmean_mc,2.),2.))*B_mc[ipt][ieta]->GetMeanError();
   Rerr_data[ipt][ieta] = sqrt(pow(2/(2-Bmean_data),2.)+pow(2*Bmean_data/pow(2-Bmean_data,2.),2.))*B_data[ipt][ieta]->GetMeanError();
   cout << ipt<<"mc error "<<Rerr_mc[ipt][ieta]<<" data error "<<Rerr_data[ipt][ieta]<<endl;
   cout << ipt<<"mc B error "<<B_mc[ipt][ieta]->GetMeanError()<<" mc B mean "<<(double)B_mc[ipt][ieta]->GetMean()<<endl;
   cout << ipt<<"data B error "<<B_data[ipt][ieta]->GetMeanError()<<" data B mean "<<(double)B_data[ipt][ieta]->GetMean()<<endl;
   corr[ieta][ipt]=R_mc[ipt][ieta]/R_data[ipt][ieta];
   correrr[ieta][ipt]=corr[ieta][ipt]*sqrt(pow(Rerr_mc[ipt][ieta]/R_mc[ipt][ieta],2.)+pow(Rerr_data[ipt][ieta]/R_data[ipt][ieta],2.));
  }
  graph[ieta] = new TGraphErrors(npt,pt[ieta],corr[ieta],pterr[ieta],correrr[ieta]);
  graphclone[ieta] = new TGraphErrors(npt,pt[ieta],corr[ieta],pterr[ieta],correrr[ieta]);
  graphclone[ieta]->Fit(fit[ieta],"R");
 }
 
 TFile *outf = new TFile(Form("Corrections/L2L3VsPtEtaBinned_alphacut_high%d_%s_etacut%d_%s.root",(int)(10*alphacut_high),algo.Data(),etacut,mode.Data()),"recreate");
 
 for(int ieta = 0; ieta < netabins; ieta++){
  graph[ieta]->SetName(Form("graph%d",ieta)); 
  graphclone[ieta]->SetName(Form("graphclone%d",ieta)); 
  graph[ieta]->Write();
  graphclone[ieta]->Write();
  fit[ieta]->Write();
  for(int ipt = 0; ipt < npt; ipt++){
   B_mc[ipt][ieta]->Write();
   B_data[ipt][ieta]->Write();
  }
 }
 outf->Close();
 
}
