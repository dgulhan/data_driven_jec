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

#include "dijeteta/utilities.h"

void extraAlpha(){
 TH1D::SetDefaultSumw2();
 TH1F::SetDefaultSumw2();
 
 TFile * file = new TFile("ntuples/photon40.root");
 TTree * ntgamma = (TTree*)file->Get("ntgammajet");
 TFile * file_mc = new TFile("ntuples/pho_pthat30.root");
 TTree * ntgamma_mc = (TTree*)file_mc->Get("ntgammajet");

 float alphalow[] = {-99,0.,0.1,0.15,0.2,0.23,0.25,0.3,0.4,0.5};
 float alphahigh[] = {0.0,0.1,0.15,0.2,0.23,0.25,0.3,0.4,0.5,1.};
 int nalpha = 9;

 // float alphalow[] = {0.1,0.15,0.2,0.23,0.25,0.3,0.4,0.5};
 // float alphahigh[] = {0.15,0.2,0.23,0.25,0.3,0.4,0.5,1.};
 // int nalpha = 7;

 TH1D * histR[nalpha];
 TH1D * histR_mc[nalpha];
 TH1D * histalpha[nalpha];
 TH1D * histalpha_mc[nalpha];
 TGraphErrors * graph;
 TGraphErrors * graph_mc;
 TGraphErrors * graph_rat;
 TGraphErrors * graph_rat_clone;
 double alpha[nalpha];
 double alphaerr[nalpha];
 double res[nalpha];
 double reserr[nalpha];
 double alpha_mc[nalpha];
 double alphaerr_mc[nalpha];
 double res_mc[nalpha];
 double reserr_mc[nalpha];
 
 double res_rat[nalpha];
 double reserr_rat[nalpha];
 
 for(int ialpha = 0; ialpha < nalpha; ialpha++){
  histR[ialpha] = new TH1D(Form("histR_%d",ialpha),"",30,0,3);
  histR_mc[ialpha] = new TH1D(Form("histR_mc_%d",ialpha),"",30,0,3);
  histalpha[ialpha] = new TH1D(Form("histalpha_%d",ialpha),"",10,alphalow[ialpha],alphahigh[ialpha]);
  histalpha_mc[ialpha] = new TH1D(Form("histalpha_mc_%d",ialpha),"",10,alphalow[ialpha],alphahigh[ialpha]);
  ntgamma->Draw(Form("res>>histR_%d",ialpha),Form("alpha>=%.2f && alpha<%.2f && phoPt1>40 && acos(cos(phoPhi1-jtphi2))>0.4 && jtpt1>10 && dphi>2.5 && abs(jteta1)<1.3 && hcalnoise",alphalow[ialpha],alphahigh[ialpha]));
  ntgamma_mc->Draw(Form("res>>histR_mc_%d",ialpha),Form("alpha>=%.2f && alpha<%.2f && phoPt1>40 && acos(cos(phoPhi1-jtphi2))>0.4 && jtpt1>10 && dphi>2.5 && abs(jteta1)<1.3",alphalow[ialpha],alphahigh[ialpha]));
  ntgamma->Draw(Form("alpha>>histalpha_%d",ialpha),Form("alpha>=%.2f && alpha<%.2f && phoPt1>40 && acos(cos(phoPhi1-jtphi2))>0.4 && jtpt1>10 && dphi>2.5 && abs(jteta1)<1.3 && hcalnoise",alphalow[ialpha],alphahigh[ialpha]));
  ntgamma_mc->Draw(Form("alpha>>histalpha_mc_%d",ialpha),Form("alpha>=%.2f && alpha<%.2f && phoPt1>40 && acos(cos(phoPhi1-jtphi2))>0.4 && jtpt1>10 && dphi>2.5 && abs(jteta1)<1.3",alphalow[ialpha],alphahigh[ialpha]));
  if(ialpha == 0){
   alpha[ialpha] = 0;
   alphaerr[ialpha] = 0;
   alpha_mc[ialpha] = 0;
   alphaerr_mc[ialpha] = 0;
  }else{
   alpha[ialpha] = histalpha[ialpha]->GetMean();
   alphaerr[ialpha] = histalpha[ialpha]->GetMeanError();
   alpha_mc[ialpha] = histalpha_mc[ialpha]->GetMean();
   alphaerr_mc[ialpha] = histalpha_mc[ialpha]->GetMeanError();
  }
  res[ialpha] = histR[ialpha]->GetMean();
  reserr[ialpha] = histR[ialpha]->GetMeanError();
  res_mc[ialpha] = histR_mc[ialpha]->GetMean();
  reserr_mc[ialpha] = histR_mc[ialpha]->GetMeanError();
  res_rat[ialpha] = res[ialpha]/res_mc[ialpha];
  reserr_rat[ialpha] = (res_mc[ialpha]/res[ialpha])*sqrt(pow(reserr[ialpha]/res[ialpha],2)+pow(reserr_mc[ialpha]/res_mc[ialpha],2));
 }
 graph = new TGraphErrors(nalpha,alpha,res,alphaerr,reserr);
 graph_mc = new TGraphErrors(nalpha,alpha_mc,res_mc,alphaerr_mc,reserr_mc);
 graph_rat = new TGraphErrors(nalpha,alpha,res_rat,alphaerr,reserr_rat);
 graph_rat_clone = new TGraphErrors(nalpha,alpha,res_rat,alphaerr,reserr_rat);
 graph->SetMarkerColor(kBlue+1);
 graph->SetLineColor(kBlue+1);
 TH1D * empty = new TH1D("empty",";p_{T,2}^{jet}/p_{T}^{#gamma};Absolute response",30,0.0001,0.4);
 empty->Fill(0.1,1000.);
 empty->SetMaximum(0.9999);
 empty->SetMinimum(0.6);
 
 empty->GetXaxis()->CenterTitle();
 empty->GetYaxis()->CenterTitle();
 TLine *line = new TLine(0.0,1.0,0.4,1.0);
 TCanvas * c1 = new TCanvas("c1","",600,600);
 empty->Draw();
 line->Draw("same");
 graph->Draw("same p");
 graph_mc->Draw("same p");
 drawText("p_{T,1}^{jet} > 10 GeV",0.5,0.9);
 drawText("p_{T}^{#gamma} > 40 GeV",0.5,0.8);
 drawText("#Delta#phi > 2.5",0.5,0.7);
 drawText("anti k_{T} R = 0.3 PF",0.22,0.25);
 
 TLegend *t8=new TLegend(0.2,0.33,0.9,0.5); 
 t8->SetFillColor(0);
 t8->SetBorderSize(0);
 t8->SetFillStyle(0);
 t8->SetTextFont(43);
 t8->SetTextSize(28); 
 t8->AddEntry(graph,"pp","p");
 t8->AddEntry(graph_mc,"PYTHIA 8","p");
 t8->Draw("same");
 
 c1->SaveAs("plots/alphadep.png");
 c1->SaveAs("plots/alphadep.pdf");
 
 TH1D * empty2 = new TH1D("empty2",";#alpha;Absolute R_{data}/R_{mc}",30,0.00001,0.4);
 empty2->SetMaximum(1.1);
 empty2->SetMinimum(0.9);
 empty2->GetXaxis()->CenterTitle();
 empty2->GetYaxis()->CenterTitle();
 TCanvas * c2 = new TCanvas("c2","",600,600);
 empty2->Draw();
 line->Draw("same");
 graph_rat->Draw("same p");
 graph_rat->SetMarkerColor(kRed+1);
 graph_rat->SetLineColor(kRed+1);
 drawText("p_{T,1}^{jet} > 10 GeV",0.5,0.9); 
 drawText("p_{T}^{#gamma} > 40 GeV",0.5,0.8);
 drawText("#Delta#phi > 2.5",0.5,0.7);
 drawText("anti k_{T} R = 0.3 PF",0.22,0.25);
 TF1 * fit = new TF1("fit","[0]*x+[1]",0.001,0.4);
 fit->SetParameters(-0.1,1.);
 graph_rat_clone->Fit(fit,"R");
 fit->SetLineColor(kBlue+1);
 fit->Draw("same");
 
 c2->SaveAs("plots/RDataOverMCVsAlpha.png");
 c2->SaveAs("plots/RDataOverMCVsAlpha.pdf");
 TFile * outf = new TFile("Corrections/extrapAlphaGammaJet_ak3PF.root","recreate");
 fit->Write();
 outf->Close();
}