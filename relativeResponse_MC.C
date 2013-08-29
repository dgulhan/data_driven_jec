#include <iostream>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom1.h"
#include "TH1F.h"
#include "TF1.h"
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

#include "plotFigure.C"
#include "akPu3PF_MC.C"

void relativeResponse_MC(){
TH2D::SetDefaultSumw2(true);
TH1D::SetDefaultSumw2();

double etacut=3;
double phicut=2.5;
char *mode="pPb";//"pPb","pp","Pbp";
char *algo="akPu3PF";
int npthat;
int pthat[10]={15,30,50,80,120,170,220,280,370,9999};
float w[9];
float w_pPb[] = {0.160697,0.0116103,0.00106034,0.000126209,2.34734e-05,3.21927e-06,6.28969e-07,2.02126e-07,5.73344e-08};//pPb
float w_Pbp[] = {0.309662,0.00815875,0.000923597,0.00010879,1.69398e-05,2.41979e-06,5.44265e-07,1.43921e-07,2.63115e-08};//pPb
float w_pp[] = {0.118313,0.00528446,0.000569802,5.90135e-05,6.35269e-06,7.11644e-07,1.21573e-07,2.65747e-08};//pPb

if(mode=="pPb" || mode=="Pbp"){
 npthat=9;
 if(mode=="pPb"){
  for(int i=0;i<npthat;i++){
   w[i]=w_pPb[i];
  }
 }
 else if(mode=="Pbp"){
  for(int i=0;i<npthat;i++){
   w[i]=w_Pbp[i];
  }
 }
}
if(mode=="pp"){
 npthat=8;
 pthat[9]=9999;
 for(int i=0; i<npthat; i++){
  w[i]=w_pp[i];
 }
}


TString infile[npthat];
HiTree *fhi[npthat];
akPu3PF_MC *t[npthat];

string dijetVars="";
string jetVars="";
dijetVars += "weight:jtpt1:jteta1:jtphi1:jtpt2:jteta2:jtphi2:jtpt3:jteta3:jtphi3:dphi"
     ":trkMax1:trkMax2"
     ":ptAverage:ptBarrel:ptProbe:etaBarrel:etaProbe:alpha:hfsum:hfp:hfm:B:pthat:refeta1:refeta2:refpt1:refpt2:refphi1:refphi2:hasforwardjet";
	 
jetVars+= "weight:jtpt:recopt:refpt:genpt:rawpt:jteta";
TFile * fnt = new TFile(Form("ntuples/ntuple_relativeResponse_MC_%s_eta%d_algo%s.root",mode,(int)etacut,algo),"recreate");

TNtuple *ntdijet = new TNtuple("ntdijet","",dijetVars.data());
TNtuple *ntjet = new TNtuple("ntjet","",jetVars.data());

TF1 *f;
if(algo=="ak3PF") f= new TF1("f","1");
if(algo=="akPu3PF"){
 f= new TF1("f","[0]/pow(x,[1])");
 f->SetParameters(0.937820,-0.0111559);
}

for(int i=0;i<npthat;i++){
 if(mode=="pPb"){
  infile[i] = Form("root://eoscms//eos/cms/store/caf/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt%d/HiForest_v77_merged01/pt%d_HP04_prod16_v77_merged_forest_0.root",pthat[i],pthat[i]);
 }
 if(mode=="Pbp"){
  if(i==0 || i==1 || i==2 || i==3 || i==4) infile[i]=Form("root://eoscms//eos/cms/store/group/phys_heavyions/dgulhan/Pbp/HP05/prod24/Hijing_Pythia_pt%d/HiForest_v84_merged02/pt%d_HP05_prod24_v84_merged_forest_0.root",pthat[i],pthat[i]);//Pbp
  if(i==5 || i==6 || i==7 || i==8) infile[i]=Form("root://eoscms//eos/cms/store/caf/user/dgulhan/Pbp/HP05/prod24/Hijing_Pythia_pt%d/HiForest_v84_merged02/pt%d_HP05_prod24_v84_merged_forest_0.root",pthat[i],pthat[i]);//Pbp
 }
 if(mode=="pp"){
  infile[i]=Form("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/inbound/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt%d/HiForest_v81_merged01/pt%d_pp2013_P01_prod22_v81_merged_forest_0.root",pthat[i],pthat[i]);
 }
 fhi[i] = new HiTree(infile[i].Data());
 t[i] = new akPu3PF_MC(infile[i].Data(),algo);
 cout<<infile[i].Data()<<endl;

 int nentries = t[i]->GetEntriesFast();
 cout<<"nentries " << nentries << " for pthat "<<pthat[i]<<endl;

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
  float ptweight = 0 ;
  float hfgenp=0;
  float hfgenm=0;
  float refeta1=-99;
  float refeta2=-99;
  float refpt1=-99;
  float refpt2=-99;
  float refphi1=-99;
  float refphi2=-99;
  bool hasForwardJet=false;
  std::vector<std::pair<float, std::pair<float, std::pair<float, std::pair<float,std::pair<float,std::pair<float,float> > > > > > > jets;
  int njet=0;
  for(int ibin = 0; ibin<npthat; ibin++){
   if(t[i]->pthat>=pthat[ibin] && t[i]->pthat<pthat[ibin+1]) ptweight = w[ibin];
  }
  for(int ijet=0;ijet<t[i]->nref;ijet++){
   if(t[i]->rawpt[ijet]<15) continue;
   if(t[i]->jtpt[ijet]<15) continue;
   if(fabs(t[i]->jteta[ijet]>etacut)) hasForwardJet=true;
   if(fabs(t[i]->jteta[ijet])>etacut) continue;
   njet++;
   float  correctedpt = (1/f->Eval(t[i]->jtpt[ijet]))*t[i]->jtpt[ijet];
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
  float  ptaverage =(ptprobe+ptbarrel)/2;
  hfsum = (fhi[i]->hiHFplusEta4 + fhi[i]->hiHFminusEta4);
  if(njet>0){
   float dijetEntry[]={ptweight,pt1,eta1,phi1,pt2,eta2,phi2,pt3,eta3,phi3,dphi,trkmax1,trkmax2,ptaverage,ptbarrel,ptprobe,etabarrel,etaprobe,(pt3/ptaverage),hfsum,fhi[i]->hiHFplusEta4,fhi[i]->hiHFminusEta4,((ptprobe-ptbarrel)/ptaverage),t[i]->pthat,refeta1,refeta2,refpt1,refpt2,refphi1,refphi2,hasForwardJet};
   ntdijet->Fill(dijetEntry);
  }
}
fhi[i]->Close();
t[i]->Close();
}

fnt->cd();
ntdijet->Write();
ntjet->Write();
fnt->Close();
cout<<"done"<<endl;
}
