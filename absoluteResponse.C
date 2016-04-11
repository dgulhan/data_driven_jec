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

#include "plotFigure.C"
#include "akPu3PF.C"

#include "L2L3ResidualWFits.h"

void absoluteResponse(TString dataset = "/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_MinBias_2/", TString infile = "MinimumBias2_HiForestAOD_101.root", TString outfname = "test.root"){ 
TH2D::SetDefaultSumw2(true);
TH1D::SetDefaultSumw2();

 int radius = 3;
 float etacut = 3;
 float phicut = 2.5;
 int docorrection = 1;
 TString mode="pp";//"pPb","pp","Pbp";
 TString algo=Form("ak%dPF",radius);
 TString infname = Form("%s/%s",dataset.Data(), infile.Data());

 L2L3Residual * L2L3Corr = new L2L3Residual(radius);

 HiTree   *fhi = new HiTree(infname.Data());
 HltTree   *fhlt = new HltTree(infname.Data());
 skimTree   *fskim = new skimTree(infname.Data());
 akPu3PF *t = new akPu3PF(infname.Data(),algo);
 EventTree *fpho = new EventTree(infname.Data());

 int nentries = t->GetEntriesFast();
 std::cout<<"nentries = "<<nentries<<std::endl;

 TString jetVars="";
 TString gammaJetVars = "";

 jetVars+= "jtpt:jteta:jtphi";
 gammaJetVars += "jtpt1:jteta1:jtphi1:jtpt2:jteta2:jtphi2:dphi"
     ":trkMax1:trkMax2"
     ":phoPt1:phoEta1:phoPhi1:phoPt2:phoEta2:phoPhi2:res:alpha:hfsum:hfp:hfm:zdcp:zdcm"
     ":hcalnoise:pileupfilter:photon40:pres40:photon30:pres30:photon50:photon60:vz:vx:vy";
 TFile * fnt = new TFile(outfname.Data(),"recreate");

 TNtuple *ntjet = new TNtuple("ntjet","",jetVars.Data());
 TNtuple *ntgammajet = new TNtuple("ntgammajet","",gammaJetVars.Data());

 for (Long64_t jentry=0; jentry<nentries;jentry++) {
  t->GetEntry(jentry); 
  fhi->GetEntry(jentry);
  fhlt->GetEntry(jentry);
  fskim->GetEntry(jentry);
  fpho->GetEntry(jentry);

  if(!(fskim->HBHENoiseFilterResultRun2Loose && fskim->pPAprimaryVertexFilter && fabs(fhi->vz)<15 && fskim->pBeamScrapingFilter)) continue;

  // if(mode=="pPb" && fhi->run>211256) continue;
  // if(mode=="Pbp" && fhi->run<=211256) continue;
  // if(fhi->run>210658) continue;


  float phopt1 = -99;
  float phopt2 = -99;
  float phopt3 = -99;
  float phoeta1 = -99;
  float phoeta2 = -99;
  float phoeta3 = -99;
  float phophi1 = -99;
  float phophi2 = -99;
  float phophi3 = -99;
    
  float pt1w=-99;
  float pt2w=-99;
  float pt3w=-99;
  float phi1w=-99;
  float phi2w=-99;
  float phi3w=-99;
  float dphiw=-99;
  float eta1w=-99;
  float eta2w=-99;
  float eta3w=-99;
  float trkmax1w=-99;
  float trkmax2w=-99;
  float ptaveragew=-99;
  
  std::vector<std::pair<float, std::pair<float, float > > > photons;
  std::vector<std::pair<float, std::pair<float, std::pair<float,float> > > > jets_corr;
  int njet = 0;
  int npho = 0;
  
  for(int ipho=0;ipho<fpho->nPho;ipho++){
   float phopt = fpho->phoEt->at(ipho);
  
   if(fabs(fpho->pho_seedTime->at(ipho)) > 3.0) continue;
   if(fabs(fpho->pho_swissCrx->at(ipho)) > 0.9) continue;
   if(!(fpho->pfcIso4->at(ipho) < 0.76  && fpho->pfnIso4->at(ipho) < (0.97 + 0.014*phopt + 0.000019*phopt*phopt) && fpho->pfpIso4->at(ipho) < (0.08 + 0.0053*phopt))) continue;
   if(fpho->phoHoverE->at(ipho)>0.05) continue;
   if(fpho->phoSigmaIEtaIEta->at(ipho) > 0.0100 ) continue;
   float phoeta = fpho->phoEta->at(ipho);
   if(fabs(phoeta)>1.44) continue;
   if(phopt < 20) continue;
   npho++;
   // corr = C_rel->GetBinContent(C_rel->FindBin(fpho->jteta[ipho]));
	// correctedpt =  correctedpt*corr;
   photons.push_back(std::make_pair(phopt,std::make_pair(phoeta, fpho->phoPhi->at(ipho))));
  }
  if(npho != 0){
   std::sort(photons.begin(),photons.end()); 
   phopt1 = photons[npho-1].first;
   phoeta1 = photons[npho-1].second.first;
   phophi1 = photons[npho-1].second.second;
   if(npho>1){  
    phopt2 = photons[npho-2].first;
    phoeta2 = photons[npho-2].second.first;
    phophi2 = photons[npho-2].second.second;
   }
  }
  
  for(int ijet=0;ijet<t->nref;ijet++){
   float jteta = t->jteta[ijet];
   float jtphi = t->jtphi[ijet];
   if(fabs(jteta)>3) continue;
   if(sqrt(pow(acos(cos(jtphi-phophi1)),2)+pow(jteta-phoeta1,2)) < 0.4) continue;
   njet++;
   float correctedpt = t->jtpt[ijet];
   // corr = C_rel->GetBinContent(C_rel->FindBin(t->jteta[ijet]));
	// correctedpt =  correctedpt*corr;
   correctedpt = L2L3Corr->get_corrected_pt(correctedpt, jteta);
   jets_corr.push_back(std::make_pair(correctedpt,std::make_pair(jteta, std::make_pair(t->jtphi[ijet],t->trackMax[ijet]))));
   ntjet->Fill(correctedpt,jteta,t->jtphi[ijet]);
  }
  // std::cout<<njet<<std::endl;
  if(njet != 0){
   std::sort(jets_corr.begin(),jets_corr.end()); 
   pt1w = jets_corr[njet-1].first;
   eta1w = jets_corr[njet-1].second.first;
   phi1w = jets_corr[njet-1].second.second.first;
   trkmax1w = jets_corr[njet-1].second.second.second;
   if(njet>1){  
    pt2w = jets_corr[njet-2].first;
    eta2w = jets_corr[njet-2].second.first;
    phi2w = jets_corr[njet-2].second.second.first;
    trkmax2w = jets_corr[njet-2].second.second.second;
	if(njet>2){
     pt3w = jets_corr[njet-3].first;
     eta3w = jets_corr[njet-3].second.first;
     phi3w = jets_corr[njet-3].second.second.first; 
    }
   }
  }
  if(fabs(eta1w)>1.3) continue;
  if(njet>=1){
   dphiw = acos(cos(phophi1-phi1w));
  }


  float gammaJetEntry[] = {pt1w,eta1w,phi1w,pt2w,eta2w,phi2w,dphiw,trkmax1w,trkmax2w,phopt1,phoeta1,phophi1,phopt2,phoeta2,phophi2,pt1w/phopt1,pt2w/phopt1,fhi->hiHFplusEta4 + fhi->hiHFminusEta4,fhi->hiHFplusEta4,fhi->hiHFminusEta4,fhi->hiZDCminus,fhi->hiZDCplus,(float)(fskim->HBHENoiseFilterResultRun2Loose),(float)(fskim->pVertexFilterCutGplus),(float)(fhlt->HLT_HISinglePhoton40_Eta3p1_v1),(float)(fhlt->HLT_HISinglePhoton40_Eta3p1_v1_Prescl),(float)(fhlt->HLT_HISinglePhoton30_Eta3p1_v1),(float)(fhlt->HLT_HISinglePhoton30_Eta3p1_v1_Prescl),(float)(fhlt->HLT_HISinglePhoton50_Eta3p1_v1),(float)(fhlt->HLT_HISinglePhoton60_Eta3p1_v1),fhi->vz,fhi->vx,fhi->vy};
  ntgammajet->Fill(gammaJetEntry);
}

ntgammajet->Write();
ntjet->Write();
fnt->Close();
}


int main(int argc, char *argv[])
{ 
  absoluteResponse(argv[1],argv[2],argv[3]);
  return 0;
}