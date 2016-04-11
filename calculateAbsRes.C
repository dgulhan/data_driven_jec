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

void calculateAbsRes(){
 TH1D::SetDefaultSumw2();
 TH1F::SetDefaultSumw2();
 
 TFile * file = new TFile("ntuples/photon40.root");
 TTree * ntgamma = (TTree*)file->Get("ntgammajet");
 TFile * file_mc[3];
 file_mc[0] = new TFile("ntuples/pho_pthat30.root");
 file_mc[1] = new TFile("ntuples/pho_pthat50.root");
 file_mc[2] = new TFile("ntuples/pho_pthat120.root");
 TTree * ntgamma_mc[3];
 for(int ifile = 0; ifile < 3; ifile++){ 
  ntgamma_mc[ifile] = (TTree*)file_mc[ifile]->Get("ntgammajet");
 }
 int index_mc_file[] = {0, 0, 1, 2};
 TFile *file_extrap = new TFile("Corrections/extrapAlphaGammaJet_ak3PF.root");
 TF1 * fit_extrap = (TF1*)file_extrap->Get("fit");
 
 int npt = 4;
 double ptlow[] = {40,45,60,140};
 double pthigh[] = {45,60,140,200};
 
 TH1D * histR[npt];
 TH1D * histR_mc[npt];
 TH1D * histpt[npt];
 TH1D * histpt_mc[npt];
 TGraphErrors * graph;
 TGraphErrors * graph_mc;
 TGraphErrors * graph_rat;
 TGraphErrors * graph_rat_clone;
 TGraphErrors * graph_rat_extrap;
 double pt[npt];
 double pterr[npt];
 double res[npt];
 double reserr[npt];
 double pt_mc[npt];
 double pterr_mc[npt];
 double res_mc[npt];
 double reserr_mc[npt];

 double res_rat[npt];
 double reserr_rat[npt];
 double res_rat_extrap[npt];
 TH1D * halpha[npt];
 for(int ipt = 0; ipt < npt; ipt++){
  histR[ipt] = new TH1D(Form("histR_%d",ipt),"",30,0,3);
  histR_mc[ipt] = new TH1D(Form("histR_mc_%d",ipt),"",30,0,3);
  histpt[ipt] = new TH1D(Form("histpt_%d",ipt),"",10,ptlow[ipt],pthigh[ipt]);
  histpt_mc[ipt] = new TH1D(Form("histpt_mc_%d",ipt),"",10,ptlow[ipt],pthigh[ipt]);
  
  ntgamma->Draw(Form("res>>histR_%d",ipt),Form("phoPt1>=%.2f && phoPt1<%.2f && phoPt1>40 && acos(cos(phoPhi1-jtphi2))>0.4 && jtpt1>10 && dphi>2.5 && abs(jteta1)<1.3 && hcalnoise && alpha < 0.2",ptlow[ipt],pthigh[ipt]));
  ntgamma_mc[index_mc_file[ipt]]->Draw(Form("res>>histR_mc_%d",ipt),Form("phoPt1>=%.2f && phoPt1<%.2f && phoPt1>40 && acos(cos(phoPhi1-jtphi2))>0.4 && jtpt1>10 && dphi>2.5 && abs(jteta1)<1.3 && alpha < 0.2 ",ptlow[ipt],pthigh[ipt]));
  ntgamma->Draw(Form("phoPt1>>histpt_%d",ipt),Form("phoPt1>=%.2f && phoPt1<%.2f && phoPt1>40 && acos(cos(phoPhi1-jtphi2))>0.4 && jtpt1>10 && dphi>2.5 && abs(jteta1)<1.3 && hcalnoise && alpha < 0.2",ptlow[ipt],pthigh[ipt]));
  ntgamma_mc[index_mc_file[ipt]]->Draw(Form("phoPt1>>histpt_mc_%d",ipt),Form("phoPt1>=%.2f && phoPt1<%.2f && phoPt1>40 && acos(cos(phoPhi1-jtphi2))>0.4 && jtpt1>10 && dphi>2.5 && abs(jteta1)<1.3 && alpha < 0.2",ptlow[ipt],pthigh[ipt]));
  pt[ipt] = histpt[ipt]->GetMean();
  pterr[ipt] = histpt[ipt]->GetMeanError();
  pt_mc[ipt] = histpt_mc[ipt]->GetMean();
  pterr_mc[ipt] = histpt_mc[ipt]->GetMeanError();

 halpha[ipt] = new TH1D(Form("halpha%d",ipt),"",30,-0.1,0.2);
 ntgamma->Draw(Form("alpha*(alpha>0)>>halpha%d",ipt),Form("phoPt1>=40 && acos(cos(phoPhi1-jtphi2))>0.4 && jtpt1>10 && dphi>2.5 && abs(jteta1)<1.3 && hcalnoise && alpha < 0.2 && phoPt1>%.2f && phoPt1<%.2f",ptlow[ipt],pthigh[ipt]));
 double meanAlpha = halpha[ipt]->GetMean();
 double extrap = fit_extrap->Eval(0.)/fit_extrap->Eval(meanAlpha);
 
  res[ipt] = histR[ipt]->GetMean();
  reserr[ipt] = histR[ipt]->GetMeanError();
  res_mc[ipt] = histR_mc[ipt]->GetMean();
  reserr_mc[ipt] = histR_mc[ipt]->GetMeanError();
  res_rat[ipt] = res[ipt]/res_mc[ipt];
  res_rat_extrap[ipt] = extrap*res_rat[ipt];
  reserr_rat[ipt] = (res_mc[ipt]/res[ipt])*sqrt(pow(reserr[ipt]/res[ipt],2)+pow(reserr_mc[ipt]/res_mc[ipt],2));
 }
 graph = new TGraphErrors(npt,pt,res,pterr,reserr);
 graph_mc = new TGraphErrors(npt,pt_mc,res_mc,pterr_mc,reserr_mc);
 graph_rat = new TGraphErrors(npt,pt,res_rat,pterr,reserr_rat);
 graph_rat_clone = new TGraphErrors(npt,pt,res_rat,pterr,reserr_rat);
 graph_rat_extrap = new TGraphErrors(npt,pt,res_rat_extrap,pterr,reserr_rat);
 graph->SetMarkerColor(kBlue+3);
 graph->SetLineColor(kBlue+3);
 TH1D * empty = new TH1D("empty",";p_{T}^{#gamma};Absolute response",30,35,180);
 empty->Fill(50,1000.);
 empty->SetMaximum(0.9999);
 empty->SetMinimum(0.7);
 
 empty->GetXaxis()->CenterTitle();
 empty->GetYaxis()->CenterTitle();
 TLine *line = new TLine(35.0,1.0,180.,1.0);
 TCanvas * c1 = new TCanvas("c1","",600,600);
 empty->Draw();
 line->Draw("same");
 graph->Draw("same p");
 graph_mc->Draw("same p");
 drawText("p_{T,1}^{jet} > 10 GeV",0.5,0.9);
 drawText("#alpha < 0.2",0.5,0.8);
 drawText("#Delta#phi > 2.5",0.5,0.7);
 drawText("anti k_{T} R = 0.3 PF",0.22,0.25);
 c1->SaveAs("plots/phoptdep.png");
 c1->SaveAs("plots/phoptdep.pdf");
 
 TH1D * empty2 = new TH1D("empty2",";p_{T}^{#gamma};Absolute R_{data}/R_{mc}",30,35,180);
 empty2->SetMaximum(1.15);
 empty2->SetMinimum(0.92);
 empty2->Fill(40,1000.);
 empty2->GetXaxis()->CenterTitle();
 empty2->GetYaxis()->CenterTitle();
 TCanvas * c2 = new TCanvas("c2","",600,600);
 empty2->Draw();
 line->Draw("same");
 graph_rat->Draw("same p");
 graph_rat->SetMarkerColor(kRed+1);
 graph_rat->SetLineColor(kRed+1);
 graph_rat_extrap->SetMarkerColor(kBlue+1);
 graph_rat_extrap->SetLineColor(kBlue+1);
 graph_rat_extrap->SetMarkerStyle(25);
 graph_rat_extrap->Draw("same p");
 drawText("p_{T,1}^{jet} > 10 GeV",0.5,0.9); 
 drawText("#alpha < 0.2",0.5,0.8);
 drawText("#Delta#phi > 2.5",0.5,0.7);
 drawText("anti k_{T} R = 0.3 PF",0.22,0.25);
 // TF1 * fit = new TF1("fit","[0]*x+[1]",0.001,0.25);
 // fit->SetParameters(-0.1,1.);
 // graph_rat_clone->Fit(fit,"R");
 // fit->SetLineColor(kBlue+1);
 // fit->Draw("same");
 
 
 TLegend *t8=new TLegend(0.2,0.7,0.5,0.9); 
 t8->SetFillColor(0);
 t8->SetBorderSize(0);
 t8->SetFillStyle(0);
 t8->SetTextFont(43);
 t8->SetTextSize(28); 
 t8->AddEntry(graph_rat,"No extrap.","p");
 t8->AddEntry(graph_rat_extrap,"With extrap.","p");
 t8->Draw("same");
 
 c2->SaveAs("plots/RDataOverMCVsPhoPt.png");
 c2->SaveAs("plots/RDataOverMCVsPhoPt.pdf");
 TFile * outf = new TFile("Corrections/absoluteResponseGammaJet_ak3PF.root","recreate");
 graph_rat_clone->SetName("response_data_over_mc");
 graph_rat_clone->Write();
 outf->Close();
 
 TH1D * empty3 = new TH1D("empty",";p_{T,1}^{jet}/p_{T}^{#gamm};Event fraction",10,0.001,1.999);
 empty3->Fill(0.,1000.);
 makePretty(empty3);
 empty3->SetMaximum(0.29999);
 empty3->SetMinimum(0.0);
 empty3->GetXaxis()->CenterTitle();
 empty3->GetYaxis()->CenterTitle();
 makePretty(empty3,2.5);
 TCanvas * c3= new TCanvas("c3","",900,900);
 makeMultiPanelCanvas(c3,2,2,0.0,0.0,0.2,0.2,0.02);
 for(int ipt = 0; ipt < npt; ipt++){ 
  c3->cd(ipt+1);
  empty3->Draw();
  histR[ipt]->Scale(1./histR[ipt]->Integral());
  histR_mc[ipt]->Scale(1./histR_mc[ipt]->Integral());
  histR_mc[ipt]->SetFillColor(18);
  histR_mc[ipt]->Draw("same hist");
  histR_mc[ipt]->SetMarkerSize(0);
  histR_mc[ipt]->Draw("same");
  histR[ipt]->SetMarkerColor(kRed+1);
  histR[ipt]->SetLineColor(kRed+1);
  histR[ipt]->Draw("same");
  drawText(Form("%d < p_{T}^{#gamma} < %d",(int)(ptlow[ipt]),(int)(pthigh[ipt])),0.25,0.9);
 }
 c3->cd(2);
 TLegend *t9=new TLegend(0.5,0.7,0.85,0.9); 
 t9->SetFillColor(0);
 t9->SetBorderSize(0);
 t9->SetFillStyle(0);
 t9->SetTextFont(43);
 t9->SetTextSize(28); 
 t9->AddEntry(histR[0],"pp 5.02 TeV","p"); 
 t9->AddEntry(histR_mc[0],"PYTHIA 8","f");
 t9->Draw("same");
 c3->SaveAs("plots/AbsResDists.png");
 c3->SaveAs("plots/AbsResDists.pdf");
}