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

void relativeResponse(){ 
TH2D::SetDefaultSumw2(true);
TH1D::SetDefaultSumw2();

double etacut=3;
double phicut=2.5;
int docorrection =1;
char *mode="pPb";//"pPb","pp","Pbp";
char *algo="akPu3PF";
TString infile;
char *dataset="jet40";//jet80,MB;
if(mode=="pPb" && dataset=="jet80") infile ="root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/pPb2013/promptReco/PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
if(mode=="pPb" && dataset=="jet40")infile ="root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/pPb2013/promptReco/PA2013_HiForest_PromptReco_JSonPPb_forestv72_HLT40_HLT60.root";
if(mode=="pPb" && dataset=="MB") infile="root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/pPb2013/promptReco/PA2013_HiForest_PromptReco_KrisztianMB_JSonPPb_forestv84.root";
// TString infile ="root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/pPb2013/promptReco/PA2013_HiForest_PromptReco_KrisztianMB_JSonPPb_forestv84.root";
if(mode=="pp" && dataset=="jet40")TString infile ="root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/pp2013/promptReco/PP2013_HiForest_PromptReco_JSon_Jet40Jet60_ppTrack_forestv84.root";
if(mode=="pp" && dataset=="jet80")TString infile ="root://eoscms//eos/cms/store/caf/user/yjlee//pp2013/promptReco/PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82.root";
if(mode=="Pbp" && dataset=="jet80")infile ="/d100/yjlee/hiForest2PPb/promptReco/Jet-Pbp-v84-JECdb/output.root";

HiTree   *fhi = new HiTree(infile.Data());
HltTree   *fhlt = new HltTree(infile.Data());
skimTree   *fskim = new skimTree(infile.Data());
akPu3PF *t = new akPu3PF(infile.Data(),algo);

TFile * fcrel3;
TH1D * C_rel;
if(docorrection){
 if(etacut==3 && mode=="pPb") fcrel3 = new TFile("Corrections/Casym_pPb_Jul17_hcalbins.root");
 if(etacut==4 && mode=="pPb")fcrel3 = new TFile("Corrections/Casym_eta4.root");
 if(etacut==3 && mode=="pp")fcrel3 = new TFile("Corrections/Casym_pp_Aug26_double_hcalbins_pt60_80.root");
 C_rel=(TH1D*)fcrel3->Get("C_asym");
}


TF1 *f;
if(algo=="ak3PF") f= new TF1("f","1");
if(algo=="akPu3PF"){
 f= new TF1("f","[0]/pow(x,[1])");
 f->SetParameters(0.937820,-0.0111559);
}

int nentries = t->GetEntriesFast();

string jetVars="";
string dijetVars = "";

jetVars+= "jtpt:jteta:jtphi";
dijetVars += "jtpt1:jteta1:jtphi1:jtpt2:jteta2:jtphi2:jtpt3:jteta3:jtphi3:dphi"
     ":trkMax1:trkMax2"
     ":ptAverage:ptBarrel:ptProbe:etaBarrel:etaProbe:alpha:hfsum:hfp:hfm:zdcp:zdcm:B"
     ":hcalnoise:pileupfilter:jet40:jet60:jet80:jet100:vz:vx:vy:hasPForwardjet:hasMForwardjet";
TFile * fnt = new TFile(Form("ntuples/ntuple_relativeResponse_%s_eta%d_corrected%d_%s_%s.root",mode,(int)etacut,docorrection,dataset,algo),"recreate");
// TFile * fnt = new TFile("ntuples/ntuple_relativeResponse_pPb_eta3_corrected_jet100_Jun03.root","recreate");

TNtuple *ntdijet = new TNtuple("ntdijet","",dijetVars.data());
TNtuple *ntjet = new TNtuple("ntjet","",jetVars.data());
TNtuple *ntdijet_corr = new TNtuple("ntdijet_corr","",dijetVars.data());

for (Long64_t jentry=0; jentry<nentries;jentry++) {
// for (Long64_t jentry=0; jentry<1000;jentry++) {
  t->GetEntry(jentry); 
  fhi->GetEntry(jentry);
  fhlt->GetEntry(jentry);
  fskim->GetEntry(jentry);

  if(!(fskim->pHBHENoiseFilter && fskim->pPAcollisionEventSelectionPA && fabs(fhi->vz)<15 && fskim->pVertexFilterCutGplus && fskim->pprimaryvertexFilter)) continue;

  if(mode=="pPb" && fhi->run>211256) continue;
  if(mode=="Pbp" && fhi->run<=211256) continue;
  // if(fhi->run>210658) continue;


  bool hasMForwardjet=false;
  bool hasPForwardjet=false;

  double pt1=-99;
  double pt2=-99;  
  double pt3=-99;
  double eta1=-99;
  double eta2=-99;  
  double eta3=-99;  
  double phi1=-99;
  double phi2=-99;
  double phi3=-99;
  double dphi=-99;
  double trkmax1 = -99;
  double trkmax2 = -99;
  double ptprobe=-99;
  double ptbarrel=-99;
  double etaprobe=-99;
  double etabarrel=-99;
  bool hasleading = false; 
  bool hassubleading = false;
  double ptaverage = -99;
  
  double pt1w=-99;
  double pt2w=-99;
  double pt3w=-99;
  double phi1w=-99;
  double phi2w=-99;
  double phi3w=-99;
  double dphiw=-99;
  double eta1w=-99;
  double eta2w=-99;
  double eta3w=-99;
  double ptprobew=-99;
  double etaprobew=-99;
  double ptbarrelw=-99;
  double etabarrelw=-99;
  double trkmax1w=-99;
  double trkmax2w=-99;
  bool hasleadingw=-99;
  bool hassubleadingw=-99;
  double ptaveragew=-99;
  
  std::vector<std::pair<double, std::pair<double, std::pair<double,double> > > > jets;
  std::vector<std::pair<double, std::pair<double, std::pair<double,double> > > > jets_corr;
  int njet = 0;
  
  for(int ijet=0;ijet<t->nref;ijet++){
   if(t->jteta[ijet]<-3) hasMForwardjet =true;
   if(t->jteta[ijet]>3) hasPForwardjet =true;
   if(fabs(t->jteta[ijet])>3) continue;
   njet++;
   double correctedpt = (f->Eval(t->jtpt[ijet]))*t->jtpt[ijet];
   jets.push_back(std::make_pair(correctedpt,std::make_pair(t->jteta[ijet], std::make_pair(t->jtphi[ijet],t->trackMax[ijet]))));
   double corr;
   if(docorrection){
    corr = C_rel->GetBinContent(C_rel->FindBin(t->jteta[ijet]));
    correctedpt = correctedpt*corr;
    jets_corr.push_back(std::make_pair(correctedpt,std::make_pair(t->jteta[ijet], std::make_pair(t->jtphi[ijet],t->trackMax[ijet]))));
   }
   ntjet->Fill(correctedpt,t->jteta[ijet],t->jtphi[ijet]);
  }
  if(njet != 0){
   std::sort(jets.begin(),jets.end());
   std::sort(jets_corr.begin(),jets_corr.end());

   pt1 = jets[njet-1].first;
   eta1 = jets[njet-1].second.first;
   phi1 = jets[njet-1].second.second.first;
   trkmax1 = jets[njet-1].second.second.second;
   pt1w = jets_corr[njet-1].first;
   eta1w = jets_corr[njet-1].second.first;
   phi1w = jets_corr[njet-1].second.second.first;
   trkmax1w = jets_corr[njet-1].second.second.second;
   if(njet>1){  
    pt2 = jets[njet-2].first;
    eta2 = jets[njet-2].second.first;
    phi2 = jets[njet-2].second.second.first;
    trkmax2 = jets[njet-2].second.second.second;
    pt2w = jets_corr[njet-2].first;
    eta2w = jets_corr[njet-2].second.first;
    phi2w = jets_corr[njet-2].second.second.first;
    trkmax2w = jets_corr[njet-2].second.second.second;
    if(njet>2){
     pt3 = jets[njet-3].first;
     eta3 = jets[njet-3].second.first;
     phi3 = jets[njet-3].second.second.first;
     pt3w = jets_corr[njet-3].first;
     eta3w = jets_corr[njet-3].second.first;
     phi3w = jets_corr[njet-3].second.second.first;
     }
   }
  }

  if(njet>=2){
   dphi = acos(cos(phi1-phi2));
   dphiw = acos(cos(phi1w-phi2w));
  }

  if ((fabs(eta1)<1.3 || fabs(eta2)<1.3) && dphi>phicut){
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
   double random = r->Rndm();
   // cout<<random<<endl;
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
  ptaverage = (ptbarrel+ptprobe)/2;

  if(docorrection){
    if ((fabs(eta1w)<1.3 || fabs(eta2w)<1.3) && dphiw>phicut){
     if(fabs(eta1w)<1.3 && fabs(eta2w)>1.3){
      etaprobew = eta2w; 
      etabarrelw = eta1w;
      ptprobew = pt2w;
      ptbarrelw = pt1w;
     }
     if(fabs(eta1w)>1.3 && fabs(eta2w)<1.3){
      etaprobew = eta1w;
      etabarrelw = eta2w;
      ptprobew = pt1w;
      ptbarrelw = pt2w;
     }
     if(fabs(eta1w)<1.3 && fabs(eta2w)<1.3){
      TRandom1 *r = new TRandom1();
      double random = r->Rndm();
      // cout<<random<<endl;
      if(random>0.5){
       etaprobew = eta1w;
       etabarrelw = eta2w;
       ptprobew = pt1w;
       ptbarrelw = pt2w;
      }else{
       etaprobew = eta2w; 
       etabarrelw = eta1w;
       ptprobew = pt2w; 
       ptbarrelw = pt1w;
      }
     }
    }
    ptaveragew =(ptprobew+ptbarrelw)/2;	
  }

  float dijetEntry[]={pt1,eta1,phi1,pt2,eta2,phi2,pt3,eta3,phi3,dphi,trkmax1,trkmax2,ptaverage,ptbarrel,ptprobe,etabarrel,etaprobe,(pt3/ptaverage),fhi->hiHFplusEta4 + fhi->hiHFminusEta4,fhi->hiHFplusEta4,fhi->hiHFminusEta4,fhi->hiZDCminus,fhi->hiZDCplus,(ptprobe-ptbarrel)/ptaverage,fskim->pHBHENoiseFilter,fskim->pVertexFilterCutGplus,fhlt->HLT_PAJet40_NoJetID_v1,fhlt->HLT_PAJet60_NoJetID_v1,fhlt->HLT_PAJet80_NoJetID_v1,fhlt->HLT_PAJet100_NoJetID_v1,fhi->vz,fhi->vx,fhi->vy,hasPForwardjet,hasMForwardjet};
  float dijetEntry_corr[] = {pt1w,eta1w,phi1w,pt2w,eta2w,phi2w,pt3w,eta3w,phi3w,dphiw,trkmax1w,trkmax2w,ptaveragew,ptbarrelw,ptprobew,etabarrelw,etaprobew,pt3w/ptaveragew,fhi->hiHFplusEta4 + fhi->hiHFminusEta4,fhi->hiHFplusEta4,fhi->hiHFminusEta4,fhi->hiZDCminus,fhi->hiZDCplus,(ptprobew-ptbarrelw)/ptaveragew,fskim->pHBHENoiseFilter,fskim->pVertexFilterCutGplus,fhlt->HLT_PAJet40_NoJetID_v1,fhlt->HLT_PAJet60_NoJetID_v1,fhlt->HLT_PAJet80_NoJetID_v1,fhlt->HLT_PAJet100_NoJetID_v1,fhi->vz,fhi->vx,fhi->vy,hasPForwardjet,hasMForwardjet};
  ntdijet->Fill(dijetEntry);
  ntdijet_corr->Fill(dijetEntry_corr);
}

ntdijet->Write();
ntjet->Write();
ntdijet_corr->Write();
fnt->Close();
}