#include <iostream>
#include <vector>
#include <algorithm>
#include <stdio.h>     
#include <math.h> 

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
#include "TCut.h"
#include "TNtuple.h"
#include "TLine.h"
#include "TF1.h"
#include "akPu3PF_MC.C"

void relativeResponse_MC(TString dataset = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt15/HiForest_v77_merged01/", TString infile = "pt15_HP04_prod16_v77_merged_forest_0.root", TString outfname = "ntuples/out_pt15_HP04_prod16_v77_merged_forest_0.root"){ 
 float etacut=3;
 float phicut=2.5;
   
 TString algo="ak3PF";
 int npthat=1; 
 // int pthat[10]={15,30,50,80,120,170,220,280,370,9999}; 

 HiTree *fhi[npthat];
 akPu3PF_MC *t[npthat];

 TString dijetVars="";
 TString jetVars="";
 TString genDijetVars="";
 dijetVars += "weight:jtpt1:jteta1:jtphi1:jtpt2:jteta2:jtphi2:jtpt3:jteta3:jtphi3:dphi"
     ":trkMax1:trkMax2"
     ":ptAverage:ptBarrel:ptProbe:etaBarrel:etaProbe:alpha:hfsum:hfp:hfm:B:pthat:refeta1:refeta2:refpt1:refpt2:refpt3:refphi1:refphi2:hasforwardjet";
 genDijetVars+="ptweight:genpt1:geneta1:genphi1:mpt1:meta1:mphi1:genpt2:geneta2:genphi2:mpt2:meta2:mphi2:pt3:eta3:phi3:mpt3:meta3:mphi3";   
	 
 jetVars+= "weight:jtpt:recopt:refpt:genpt:rawpt:jteta";
 TFile * fnt = new TFile(outfname.Data(),"recreate");

 TNtuple *ntdijet = new TNtuple("ntdijet","",dijetVars.Data());
 TNtuple *ntgendijet = new TNtuple("ntgendijet","",genDijetVars.Data());
 TNtuple *ntjet = new TNtuple("ntjet","",jetVars.Data());

 for(int i=0;i<1;i++){
 
  TString infname = Form("%s/%s",dataset.Data(), infile.Data());

  fhi[i] = new HiTree(infname.Data());
  t[i] = new akPu3PF_MC(infname.Data(),algo);
  std::cout<<infname.Data()<<std::endl;

  int nentries = t[i]->GetEntriesFast();
  // std::cout<<"nentries " << nentries << " for pthat "<<pthat[i]<<std::endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  // for (Long64_t jentry=0; jentry<100;jentry++) {
   t[i]->GetEntry(jentry); 
   fhi[i]->GetEntry(jentry);
   float pt1=-99;
   float pt2=-99;  
   float pt3=-99;
   float eta1=-99;
   float eta2=-99;  
   float eta3=-99;  
   float phi1=-99;
   float phi2=-99;
   float phi3=-99;
   float dphi=-99;
   float trkmax1=-99;
   float trkmax2=-99;
   float ptprobe=-99;
   float ptbarrel=-99;
   float etaprobe=-99;
   float etabarrel=-99;
   float hfsum = -99;
   float c = 1;
   float ptweight = 1 ;
   float hfgenp=0;
   float hfgenm=0;
   float refeta1=-99;
   float refeta2=-99;
   float refpt1=-99;
   float refpt2=-99;
   float refphi1=-99;
   float refphi2=-99;
   float refpt3=-99;
  
  
   float genpt1=-99;
   float genpt2=-99;  
   float genpt3=-99;
   float geneta1=-99;
   float geneta2=-99;  
   float geneta3=-99;  
   float genphi1=-99;
   float genphi2=-99;
   float genphi3=-99;
   float mpt1=-99;
   float mpt2=-99;  
   float mpt3=-99;
   float meta1=-99;
   float meta2=-99;  
   float meta3=-99;  
   float mphi1=-99;
   float mphi2=-99;
   float mphi3=-99;
   bool hasForwardJet=false;
   std::vector<std::pair<float, std::pair<float, std::pair<float, std::pair<float,std::pair<float,std::pair<float,float> > > > > > > jets;
   std::vector<std::pair<float, std::pair<float, std::pair<float, std::pair<float,std::pair<float,float > > > > > > genjets;
   int njet=0;

   for(int ijet=0;ijet<t[i]->nref;ijet++){
    if(fabs(t[i]->jteta[ijet]>etacut)) hasForwardJet=true;
    if(fabs(t[i]->jteta[ijet])>etacut) continue;
    njet++;
    float  correctedpt = t[i]->jtpt[ijet];
    // double correctedpt = (f->Eval(t[i]->jtpt[ijet]))*t[i]->jtpt[ijet];
    jets.push_back(std::make_pair(correctedpt,std::make_pair(t[i]->jteta[ijet], std::make_pair(t[i]->jtphi[ijet],std::make_pair(t[i]->trackMax[ijet],std::make_pair(t[i]->refeta[ijet],std::make_pair(t[i]->refpt[ijet],t[i]->refphi[ijet])))))));
    ntjet->Fill(ptweight,correctedpt,t[i]->jtpt[ijet],t[i]->refpt[ijet],t[i]->genpt[ijet],t[i]->jteta[ijet]);
   }
   if(njet>0){
    std::sort(jets.begin(),jets.end());
    pt1 = jets[njet-1].first;
    eta1 = jets[njet-1].second.first;
    phi1 = jets[njet-1].second.second.first;
    trkmax1 = jets[njet-1].second.second.second.first;
    refeta1 = jets[njet-1].second.second.second.second.first;
    refpt1 = jets[njet-1].second.second.second.second.second.first;
    refphi1 = jets[njet-1].second.second.second.second.second.second;
    if(njet>1){
     pt2 = jets[njet-2].first;
     eta2 = jets[njet-2].second.first;
     phi2 = jets[njet-2].second.second.first;
     trkmax2 = jets[njet-2].second.second.second.first;
     refeta2 = jets[njet-2].second.second.second.second.first;
     refpt2 = jets[njet-2].second.second.second.second.second.first;
     refphi2 = jets[njet-2].second.second.second.second.second.second;
     if(njet>2){
      pt3 = jets[njet-3].first;
      eta3 = jets[njet-3].second.first;
      phi3 = jets[njet-3].second.second.first;
      refpt3 = jets[njet-3].second.second.second.second.second.first;
     }
    }
    }
    if(njet>=2) dphi = acos(cos(phi1-phi2));
    if (fabs(eta1)<1.3 || fabs(eta2)<1.3){ 
    if(fabs(eta1)<1.3 && fabs(eta2)>1.3){
     etaprobe = eta2;
     etabarrel = eta1;
     ptprobe = pt2;
     ptbarrel = pt1;
    }
    if(fabs(eta1)>1.3 && fabs(eta2)<1.3){
     etaprobe = eta1;
     etabarrel = eta2;
     ptprobe = pt1;
     ptbarrel = pt2;
    }
    if(fabs(eta1)<1.3 && fabs(eta2)<1.3){
     TRandom1 *r = new TRandom1();
     float  random = r->Rndm();
     if(random>0.5){
      etaprobe = eta1;
      etabarrel = eta2;
      ptprobe = pt1;
      ptbarrel = pt2;
     }else{
      etaprobe = eta2; 
      etabarrel = eta1;
      ptprobe = pt2; 
      ptbarrel = pt1;
     }
    }
   }
  
   int ngenjet = 0;
   for(int ijet=0;ijet<t[i]->ngen;ijet++){
    if(fabs(t[i]->geneta[ijet])>etacut) continue;
    float pt = t[i]->genpt[ijet];
    float eta = t[i]->geneta[ijet];
    float phi = t[i]->genphi[ijet];
    float mpt = -99;
    float meta = -99;
    float mphi = -99;
    for(int ireco = 0; ireco < njet; ireco++){
     if(fabs(jets[ireco].first/pt-1)<0.9 && fabs(jets[ireco].second.first - eta)<0.4 && fabs(jets[ireco].second.second.first - phi)<0.4){
	  mpt = jets[ireco].first;
	  meta = jets[ireco].second.first;
	  mphi = jets[ireco].second.second.first;
	 }
    }
    ngenjet++;
    genjets.push_back(std::make_pair(pt,std::make_pair(eta, std::make_pair(phi,std::make_pair(mpt,std::make_pair(meta,mphi))))));
   }
   if(ngenjet>0){
    std::sort(genjets.begin(),genjets.end());
    genpt1 = genjets[ngenjet-1].first;
    geneta1 = genjets[ngenjet-1].second.first;
    genphi1 = genjets[ngenjet-1].second.second.first;
    mpt1 = genjets[ngenjet-1].second.second.second.first;
    meta1 = genjets[ngenjet-1].second.second.second.second.first;
    mphi1 = genjets[ngenjet-1].second.second.second.second.second;
    if(ngenjet>1){
     genpt2 = genjets[ngenjet-2].first;
     geneta2 = genjets[ngenjet-2].second.first;
     genphi2 = genjets[ngenjet-2].second.second.first;
     mpt2 = genjets[ngenjet-2].second.second.second.first;
     meta2 = genjets[ngenjet-2].second.second.second.second.first;
     mphi2 = genjets[ngenjet-2].second.second.second.second.second;
     if(ngenjet>2){
      genpt3 = genjets[ngenjet-3].first;
      geneta3 = genjets[ngenjet-3].second.first;
      genphi3 = genjets[ngenjet-3].second.second.first;
      mpt3 = genjets[ngenjet-3].second.second.second.first;
      meta3 = genjets[ngenjet-3].second.second.second.second.first;
      mphi3 = genjets[ngenjet-3].second.second.second.second.second;
     }
    }
   }
   float  ptaverage =(ptprobe+ptbarrel)/2;
   hfsum = (fhi[i]->hiHFplusEta4 + fhi[i]->hiHFminusEta4);
   if(njet>0){
    float dijetEntry[]={ptweight,pt1,eta1,phi1,pt2,eta2,phi2,pt3,eta3,phi3,dphi,trkmax1,trkmax2,ptaverage,ptbarrel,ptprobe,etabarrel,etaprobe,(pt3/ptaverage),hfsum,fhi[i]->hiHFplusEta4,fhi[i]->hiHFminusEta4,((ptprobe-ptbarrel)/ptaverage),t[i]->pthat,refeta1,refeta2,refpt1,refpt2,refpt3,refphi1,refphi2,((float)hasForwardJet)};
    ntdijet->Fill(dijetEntry);
   }
   if(ngenjet>0){
    float genDijetEntry[]={ptweight,genpt1,geneta1,genphi1,mpt1,meta1,mphi1,genpt2,geneta2,genphi2,mpt2,meta2,mphi2,genpt3,geneta3,genphi3,mpt3,meta3,mphi3};   
    ntgendijet->Fill(genDijetEntry);
   }
  }
  fhi[i]->Close();
  t[i]->Close();
 }

 fnt->cd(); 
 ntgendijet->Write();
 ntdijet->Write();
 ntjet->Write();
 fnt->Close();
 std::cout<<"done"<<std::endl;
}

int main(int argc, char *argv[])
{ 
  relativeResponse_MC(argv[1],argv[2],argv[3]);
  return 0;
}