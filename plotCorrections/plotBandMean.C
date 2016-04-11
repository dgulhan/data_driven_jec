#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TLine.h"

#include "../dijeteta/utilities.h"

void plotBandMean(){
 TH1D::SetDefaultSumw2();
int etacut = 4;
 
 TFile *file = new TFile(Form("../Corrections/L2L3VsPtEtaBinned_alphacut_high2_ak3PF_etacut%d.root",etacut));
 cout << Form("../Corrections/L2L3VsPtEtaBinned_alphacut_high2_ak3PF_etacut%d.root",etacut)<<endl;
 int npt = 7;
int netabins=64; 
// int netabins=1;
double etabins_double_hcalbins[netabins];
 
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
					  
double etabins_hcalbins4[]= {-4,      -3.664,  -3.314,  -2.964, -2.853,
                     -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830,
                     -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
                     -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609,
                     -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,  0.000,
                      0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,
                      0.696,  0.783,  0.879,  0.957,  1.044,  1.131,  1.218,
                      1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,
                      1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  2.853,
                      2.964,  3.314,  3.664, 4};
 double etabins_double_hcalbins4[]={-4,      -3.664,  -3.314,  -2.964, 
                  -2.650,  -2.322,  -2.043,  -1.830,
                  -1.653,  -1.479, -1.305, 
                  -1.131,  -0.957, -0.783, -0.609,
                   -0.435, -0.261,  -0.087,
                   0.087,  0.261,  0.435,  0.609,
                   0.783,  0.957,  1.131, 
                   1.305,  1.479,  1.653,  1.830,
                   2.043,  2.322,  2.650, 
                   2.964,  3.314,  3.664, 4};
				   
double etabins_double_hcalbins3[]= {-3,-2.500,  -2.172,-1.740, -1.392, -1.044,  -0.696,  -0.348,  0.000,  0.348,   0.696,   1.044,1.392,   1.740,   2.172,   2.500,   3};
if(etacut==3){
 // netabins=16;
 netabins=58;
 for(int ieta = 0; ieta<netabins+1; ieta++){
  // etabins_double_hcalbins[ieta] = etabins_double_hcalbins3[ieta];
  etabins_double_hcalbins[ieta] = etabins_hcalbins[ieta];
 }
}
if(etacut==4){				
 for(int ieta = 0; ieta<netabins+1; ieta++){
  // etabins_double_hcalbins[ieta] = etabins_double_hcalbins4[ieta];
  etabins_double_hcalbins[ieta] = etabins_hcalbins4[ieta];
 }  
}
 
 int ptlow[]= {30,60,80,95,120,150,200}; 
 int pthigh[]={60,80,95,120,150,200,300};
             
 TH1D *B_mc[npt][netabins];
 TH1D *B_data[npt][netabins];
 
 TLegend *t8=new TLegend(0.05,0.25,0.9,0.4); 
 t8->SetFillColor(0);
 t8->SetBorderSize(0);
 t8->SetFillStyle(0);
 t8->SetTextFont(43);
 t8->SetTextSize(28); 
 
 for(int ipt = 0; ipt < npt; ipt++){
  for(int ieta = 0; ieta < netabins; ieta++){
   B_mc[ipt][ieta] = (TH1D*)file->Get(Form("B_mc_%d_%d",ipt,ieta));
   B_data[ipt][ieta] = (TH1D*)file->Get(Form("B_data_%d_%d",ipt,ieta));
   B_mc[ipt][ieta]->SetFillColor(18);
   // B_mc[ipt][ieta]->SetLineColor(18); 
   B_mc[ipt][ieta]->SetMarkerSize(0);
   B_data[ipt][ieta]->SetMarkerColor(kRed+1);
   B_data[ipt][ieta]->SetLineColor(kRed+1);
   B_mc[ipt][ieta]->Scale(1./B_mc[ipt][ieta]->Integral());
   B_data[ipt][ieta]->Scale(1./B_data[ipt][ieta]->Integral());
  }
 }
 t8->AddEntry(B_mc[0][0],"PYTHIA 8","f");
 t8->AddEntry(B_data[0][0],"pp","p");
 
 TH1D * empty = new TH1D("empty",";B;Event fraction",10,-1.5,1.5);
 empty->Fill(0.,1000.);
 makePretty(empty);
 empty->SetMaximum(0.2);
 empty->SetMinimum(0.0);
 TCanvas * c1[netabins];
 for(int ieta = 0; ieta < netabins; ieta++){
  c1[ieta] = new TCanvas(Form("c1_%d",ieta),"",1200,800);
  makeMultiPanelCanvas(c1[ieta],4,2,0.0,0.0,0.2,0.2,0.02);
  c1[ieta]->cd(1);
  drawText(Form("%.3f < #eta_{probe} < %.3f",etabins_double_hcalbins[ieta],etabins_double_hcalbins[ieta+1]),0.05,0.5,28);
  t8->Draw("same");
  for(int ipt = 0; ipt<npt; ipt++){
   c1[ieta]->cd(ipt+2);
   empty->Draw();
   B_mc[ipt][ieta]->Draw("same"); 
   B_mc[ipt][ieta]->Draw("same hist");
   B_data[ipt][ieta]->Draw("same");
   if(ipt == 3) drawText(Form("%d < p_{T,ave} < %d", ptlow[ipt],pthigh[ipt]),0.25,0.9,28);
   else drawText(Form("%d < p_{T,ave} < %d", ptlow[ipt],pthigh[ipt]),0.05,0.9,28);
   c1[ieta]->cd(ipt+2)->RedrawAxis();
  }
  c1[ieta]->SaveAs(Form("balance%d.pdf",ieta));
  c1[ieta]->SaveAs(Form("balance%d.png",ieta));
 }
 
 TCanvas * c2 = new TCanvas("c2","",1600,1600);
 makeMultiPanelCanvas(c2,4,4,0.0,0.0,0.25,0.25,0.02);
 TGraphErrors * graph[netabins];
 TF1 * fit[netabins];
 TH1D *empty2 = new TH1D("empty2",";p_{T,ave};R_{mc}/R_{data}",50,20.0001,249.999);
 empty2->Fill(100,1000.);
 if(etacut==4){
 empty2->SetMaximum(1.2999);
 empty2->SetMinimum(0.50001);
 }
 if(etacut==3){
 
 empty2->SetMaximum(1.0999);
 empty2->SetMinimum(0.90001);
 }
 makePretty(empty2,5.0);
 TLine *line = new TLine(20,1.,250,1.);
 line->SetLineStyle(2);
 for(int ieta = 0; ieta < 16; ieta++){
  c2->cd(ieta+1);
  empty2->Draw();
  graph[ieta] = (TGraphErrors*)file->Get(Form("graph%d",ieta));
  fit[ieta] = (TF1*)file->Get(Form("fit%d",ieta));
  graph[ieta]->Draw("same p");
  fit[ieta]->Draw("same");
  line->Draw("same");
  if(ieta == 0 || ieta == 4 || ieta == 8 || ieta == 12)drawText(Form("%.3f < #eta_{probe} < %.3f",etabins_double_hcalbins[ieta],etabins_double_hcalbins[ieta+1]),0.25,0.9,28);
  else drawText(Form("%.3f < #eta_{probe} < %.3f",etabins_double_hcalbins[ieta],etabins_double_hcalbins[ieta+1]),0.05,0.9,28);
 }
 c2->SaveAs(Form("averagesAndFits%d.png",etacut)); 
 c2->SaveAs(Form("averagesAndFits%d.pdf",etacut)); 
 // if(etacut==4){
 TCanvas * c3 = new TCanvas("c3","",1600,1600);
 makeMultiPanelCanvas(c3,4,4,0.0,0.0,0.25,0.25,0.02);
 TGraphErrors * graph3[netabins];
 TF1 * fit3[netabins];
 for(int ieta = 16; ieta < 32; ieta++){
  c3->cd(ieta-16+1);
  empty2->Draw();
  graph3[ieta] = (TGraphErrors*)file->Get(Form("graph%d",ieta));
  fit3[ieta] = (TF1*)file->Get(Form("fit%d",ieta));
  graph3[ieta]->Draw("same p");
  fit3[ieta]->Draw("same");
  line->Draw("same");
  drawText(Form("%.3f < #eta_{probe} < %.3f",etabins_double_hcalbins[ieta],etabins_double_hcalbins[ieta+1]),0.25,0.9,28);
  }
 c3->SaveAs(Form("averagesAndFits%d_panels2.png",etacut)); 
 c3->SaveAs(Form("averagesAndFits%d_panels2.pdf",etacut)); 
 // } 
 TCanvas * c4 = new TCanvas("c4","",1600,2000);
 makeMultiPanelCanvas(c4,4,4,0.0,0.0,0.25,0.25,0.02);
 TGraphErrors * graph4[netabins];
 TF1 * fit4[netabins];
 for(int ieta = 32; ieta < 48; ieta++){
  c4->cd(ieta-32+1);
  empty2->Draw();
  graph4[ieta] = (TGraphErrors*)file->Get(Form("graph%d",ieta));
  fit4[ieta] = (TF1*)file->Get(Form("fit%d",ieta));
  graph4[ieta]->Draw("same p");
  fit4[ieta]->Draw("same");
  line->Draw("same");
  drawText(Form("%.3f < #eta_{probe} < %.3f",etabins_double_hcalbins[ieta],etabins_double_hcalbins[ieta+1]),0.25,0.9,28);
  }
 c4->SaveAs(Form("averagesAndFits%d_panels3.png",etacut)); 
 c4->SaveAs(Form("averagesAndFits%d_panels3.pdf",etacut)); 
 if(etacut==4){
 TCanvas * c5 = new TCanvas("c5","",1600,2000);
 makeMultiPanelCanvas(c5,4,5,0.0,0.0,0.25,0.25,0.02);
 TGraphErrors * graph5[netabins];
 TF1 * fit5[netabins];
 for(int ieta = 48; ieta < netabins; ieta++){
  c5->cd(ieta-48+1);
  empty2->Draw();
  graph4[ieta] = (TGraphErrors*)file->Get(Form("graph%d",ieta));
  fit4[ieta] = (TF1*)file->Get(Form("fit%d",ieta));
  graph4[ieta]->Draw("same p");
  fit4[ieta]->Draw("same");
  line->Draw("same");
  drawText(Form("%.3f < #eta_{probe} < %.3f",etabins_double_hcalbins[ieta],etabins_double_hcalbins[ieta+1]),0.25,0.9,28);
  }
 c5->SaveAs(Form("averagesAndFits%d_panels4.png",etacut)); 
 c5->SaveAs(Form("averagesAndFits%d_panels4.pdf",etacut)); 
 }
}