#include <iostream>
#include <vector>
#include <algorithm>

#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom3.h"
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

void relativeResponse(TString dataset = "/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_MinBias_2/", TString infile = "MinimumBias2_HiForestAOD_101.root", TString outfname = "test.root"){ 
TH2D::SetDefaultSumw2(true);
TH1D::SetDefaultSumw2();

int radius = 3; 
float etacut = 3;
float cutband = 0.0;
float etacutcorr = etacut;
bool doShift = true;
float shift = -0.465;
if(doShift){
 etacutcorr = 4;
 shift = -0.465;
}
float phicut = 2.5;
int docorrection = 1;
TString mode="pp";//"pPb","pp","Pbp";
TString algo=Form("ak%dPF",radius);
TString infname = Form("%s/%s",dataset.Data(), infile.Data());
L2L3Residual * L2L3Corr = new L2L3Residual(radius,((int)etacutcorr));
HiTree   *fhi = new HiTree(infname.Data()); 
HltTree   *fhlt = new HltTree(infname.Data());
skimTree   *fskim = new skimTree(infname.Data());
akPu3PF *t = new akPu3PF(infname.Data(),algo);

TFile * fcrel3;
TH1D * C_rel;

TFile *fileResJESSys[3];
TH1D *histJESVar[3];
TString corrstring[] = {"pt55_75_lowerpt","pt75_100_jet60","pt100_400_jet80"};
if(etacutcorr==4){
 for(int ifile = 0; ifile < 3; ifile++){
  fileResJESSys[ifile] = new TFile(Form("Casym_pp_hcalbins4_algo_ak3PF_%s_corrv5_eta4_alphahigh_20_phicut250_etacut4.root",corrstring[ifile].Data()));
  histJESVar[ifile] = (TH1D*)fileResJESSys[ifile]->Get("R_data_MC_corr");
 }
}
if(etacutcorr==3){
 for(int ifile = 0; ifile < 3; ifile++){
  fileResJESSys[ifile] = new TFile(Form("Casym_pp_hcalbins_algo_ak3PF_%s_corrv5_eta3_alphahigh_20_phicut250_etacut0.root",corrstring[ifile].Data()));
  histJESVar[ifile] = (TH1D*)fileResJESSys[ifile]->Get("R_data_MC_corr");
 }
}

if(0){
 if(etacut==3 && mode=="pPb" && algo=="akPu3PF") fcrel3 = new TFile("Corrections/Casym_pPb_double_hcalbins_algo_akPu3PF_pt100_140_jet80_alphahigh_20_phicut250.root");
 if(etacut==3 && mode=="pPb" && algo=="ak3PF") fcrel3 = new TFile("Corrections/Casym_pPb_double_hcalbins_algo_ak3PF_pt100_140_jet80_alphahigh_20_phicut250.root");
 if(etacut==4 && mode=="pPb")fcrel3 = new TFile("Corrections/Casym_eta4.root");
 if(etacut==3 && mode=="pp" && algo=="ak3PF")fcrel3 = new TFile("Casym_pp_double_hcalbins_algo_ak3PF_pt100_140_jet80_alphahigh_20_phicut250.root");
 if(etacut==3 && mode=="pp" && algo=="akPu3PF")fcrel3 = new TFile("Corrections/Casym_pp_double_hcalbins_algo_akPu3PF_pt100_140_jet80_alphahigh_20_phicut250.root");
 if(etacut==3 && mode=="pp" && algo=="AK4PF") fcrel3 = new TFile("Casym_pp_double_hcalbins_algo_AK4PF_pt100_140_jet80_alphahigh_20_phicut250.root");
 if(etacut==3 && mode=="pp" && algo=="ak4PF") fcrel3 = new TFile("Casym_pp_double_hcalbins_algo_ak4PF_pt100_140_jet80_alphahigh_20_phicut250.root");
 if(etacut==3 && mode=="pp" && algo=="ak5PF") fcrel3 = new TFile("Casym_pp_double_hcalbins_algo_ak5PF_pt100_140_jet80_alphahigh_20_phicut250.root");
 C_rel=(TH1D*)fcrel3->Get("C_asym");
}


TF1 *f;
if(algo=="ak3PF") f= new TF1("f","1");
if(algo=="akPu3PF"){
 if(mode=="pp") f= new TF1("f","1");
 if(mode=="pPb" || mode=="Pbp"){
  f= new TF1("f","[0]/pow(x,[1])");
  f->SetParameters(0.937820,-0.0111559);
 }
}

int nentries = t->GetEntriesFast();
std::cout<<"nentries = "<<nentries<<std::endl;

TString jetVars="";
TString dijetVars = "";
TString debugVars = "";

debugVars += "jtpt1:jteta1:trkMax1:trkSum1:chargedN1:rawpt1:eSum1";
jetVars+= "jtpt:jteta:jtphi";
dijetVars += "jtpt1:jteta1:jtphi1:jtpt2:jteta2:jtphi2:jtpt3:jteta3:jtphi3:dphi"
     ":trkMax1:trkMax2"
     ":ptAverage:ptBarrel:ptProbe:etaBarrel:etaProbe:alpha:hfsum:hfp:hfm:zdcp:zdcm:B"
     ":hcalnoise:pileupfilter:jet40:presjet40:jet60:presjet60:jet80:jet100:vz:vx:vy:hasPForwardjet:hasMForwardjet:random:jtpt1up:jtpt1down:jtpt2up:jtpt2down";
TFile * fnt = new TFile(outfname.Data(),"recreate");

// TFile * fnt = new TFile("ntuples/ntuple_relativeResponse_pPb_eta3_corrected_jet100_Jun03.root","recreate");

TNtuple *ntdijet = new TNtuple("ntdijet","",dijetVars.Data());
TNtuple *ntjet = new TNtuple("ntjet","",jetVars.Data());
TNtuple *ntdijet_corr = new TNtuple("ntdijet_corr","",dijetVars.Data());
TNtuple *ntdebug = new TNtuple("ntdebug","",debugVars.Data());
TRandom3 *r = new TRandom3();

for (Long64_t jentry=0; jentry<nentries;jentry++) {
// for (Long64_t jentry=0; jentry<1000;jentry++) {
  t->GetEntry(jentry); 
  fhi->GetEntry(jentry);
  fhlt->GetEntry(jentry);
  fskim->GetEntry(jentry);
  if(!(fskim->HBHENoiseFilterResultRun2Loose && fskim->pPAprimaryVertexFilter && fabs(fhi->vz)<15 && fskim->pBeamScrapingFilter)) continue;
  if(mode=="pPb" && fhi->run>211256) continue;
  if(mode=="Pbp" && fhi->run<=211256) continue;
  // if(fhi->run>210658) continue;


  bool hasMForwardjet=false;
  bool hasPForwardjet=false;

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
  float trkmax1 = -99;
  float trkSum1 = -99;
  float chargedN1 = -99;
  float rawpt1 = -99;
  float eSum1 = -99;
  float trkmax2 = -99;
  float ptprobe=-99;
  float ptbarrel=-99;
  float etaprobe=-99;
  float etabarrel=-99;
  bool hasleading = false; 
  bool hassubleading = false;
  float ptaverage = -99;
  
  float pt1wup=-99;
  float pt2wup=-99;
  float pt1wdown=-99;
  float pt2wdown=-99;
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
  float ptprobew=-99;
  float etaprobew=-99;
  float ptbarrelw=-99;
  float etabarrelw=-99;
  float trkmax1w=-99;
  float trkmax2w=-99;
  bool hasleadingw=-99;
  bool hassubleadingw=-99;
  float ptaveragew=-99;
  
  std::vector<std::pair<float, std::pair<float, std::pair<float,std::pair<float,std::pair<float,std::pair<float,std::pair<float,float> > > > > > > > jets;
  std::vector<std::pair<float, std::pair<float, std::pair<float,float> > > > jets_corr;
  std::vector<std::pair<float, std::pair<float, std::pair<float,float> > > > rawjets;
  int njet = 0;
  
  for(int ijet=0;ijet<t->nref;ijet++){
   float jteta = t->jteta[ijet];
   if(jteta<-3) hasMForwardjet =true;
   if(jteta>3) hasPForwardjet =true;
   if(fabs(jteta+shift)>etacut-cutband) continue;
   // if(t->rawpt[ijet]<15) continue;

   njet++;
   // std::cout<<"event with jet"<<njet<<std::endl;
   float correctedpt = t->jtpt[ijet];
   jets.push_back(std::make_pair(correctedpt,std::make_pair(jteta, std::make_pair(t->jtphi[ijet], std::make_pair(t->trackMax[ijet], std::make_pair(t->trackSum[ijet], std::make_pair(t->chargedN[ijet], std::make_pair(t->rawpt[ijet],t->eSum[ijet]))))))));
   rawjets.push_back(std::make_pair(t->rawpt[ijet],std::make_pair(jteta, std::make_pair(t->jtphi[ijet],t->trackMax[ijet]))));

   float corr;
   if(docorrection){
    correctedpt = L2L3Corr->get_corrected_pt(correctedpt, jteta);
    jets_corr.push_back(std::make_pair(correctedpt,std::make_pair(jteta, std::make_pair(t->jtphi[ijet],t->trackMax[ijet]))));
   }
   ntjet->Fill(correctedpt,jteta,t->jtphi[ijet]);
  }
  if(njet != 0){
   std::sort(jets.begin(),jets.end()); 
   std::sort(jets_corr.begin(),jets_corr.end()); 

   pt1 = jets[njet-1].first;
   eta1 = jets[njet-1].second.first+shift;
   phi1 = jets[njet-1].second.second.first;
   trkmax1 = jets[njet-1].second.second.second.first;
   trkSum1 = jets[njet-1].second.second.second.second.first;
   chargedN1 = jets[njet-1].second.second.second.second.second.first;
   rawpt1 = jets[njet-1].second.second.second.second.second.second.first;
   eSum1 = jets[njet-1].second.second.second.second.second.second.second;
   if(docorrection){
    pt1w = jets_corr[njet-1].first;
    eta1w = jets_corr[njet-1].second.first+shift;
	int ifile = 0;
	if(pt1w < 85) ifile = 0;
	if(pt1w > 85 && pt1w < 100) ifile = 1;
	if(pt1w > 100 && pt1w < 400) ifile = 2;
	pt1wup = (sqrt(pow(fabs(histJESVar[ifile]->GetBinContent(histJESVar[ifile]->FindBin(eta1w-shift))-1),2)+0.02*0.02)+1)*pt1w;
	pt1wdown = (1-sqrt(pow(fabs(histJESVar[ifile]->GetBinContent(histJESVar[ifile]->FindBin(eta1w-shift))-1),2)+0.02*0.02))*pt1w;
    phi1w = jets_corr[njet-1].second.second.first;
    trkmax1w = jets_corr[njet-1].second.second.second;
   }
   if(njet>1){  
    pt2 = jets[njet-2].first;
    eta2 = jets[njet-2].second.first+shift;
    phi2 = jets[njet-2].second.second.first;
    trkmax2 = jets[njet-2].second.second.second.first;
	if(docorrection){
    pt2w = jets_corr[njet-2].first;
    eta2w = jets_corr[njet-2].second.first+shift;
	int ifile = 0;
	if(pt2w < 85) ifile = 0;
	if(pt2w > 85 && pt2w < 100) ifile = 1;
	if(pt2w > 100 && pt2w < 400) ifile = 2;
	pt2wup = (sqrt(pow(fabs(histJESVar[ifile]->GetBinContent(histJESVar[ifile]->FindBin(eta2w-shift))-1),2)+0.02*0.02)+1)*pt2w;
	pt2wdown = (1-sqrt(pow(fabs(histJESVar[ifile]->GetBinContent(histJESVar[ifile]->FindBin(eta2w-shift))-1),2)+0.02*0.02))*pt2w;
    phi2w = jets_corr[njet-2].second.second.first;
    trkmax2w = jets_corr[njet-2].second.second.second;
	}
    if(njet>2){
     pt3 = jets[njet-3].first;
     eta3 = jets[njet-3].second.first+shift;
     phi3 = jets[njet-3].second.second.first;
	 if(docorrection){
     pt3w = jets_corr[njet-3].first;
     eta3w = jets_corr[njet-3].second.first+shift;
     phi3w = jets_corr[njet-3].second.second.first;
	 }
     }
   }
  }
  
  if(njet>=2){
   dphi = acos(cos(phi1-phi2));
   dphiw = acos(cos(phi1w-phi2w));
  }
  
  float randnum = -99;
  if ((fabs(eta1)<1.3 || fabs(eta2)<1.3) && dphi>phicut){
  if(fabs(eta1)<1.3 && fabs(eta2)>=1.3){
   etaprobe = eta2; 
   etabarrel = eta1;
   ptprobe = pt2;
   ptbarrel = pt1;
  }
  if(fabs(eta1)>=1.3 && fabs(eta2)<1.3){
   etaprobe = eta1;
   etabarrel = eta2;
   ptprobe = pt1;
   ptbarrel = pt2;
  }
  if(fabs(eta1)<1.3 && fabs(eta2)<1.3){
   float random = r->Rndm();
   randnum = random;
   // std::cout<<random<<std::endl;
   // if(random>0.47){
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
  float randnumw = -99;

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
      float random = r->Rndm();
	  randnumw = random;
      // std::cout<<random<<std::endl;
      // if(random>0.47){
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

  float dijetEntry[]={pt1,eta1,phi1,pt2,eta2,phi2,pt3,eta3,phi3,dphi,trkmax1,trkmax2,ptaverage,ptbarrel,ptprobe,etabarrel,etaprobe,(pt3/ptaverage),fhi->hiHFplusEta4 + fhi->hiHFminusEta4,fhi->hiHFplusEta4,fhi->hiHFminusEta4,fhi->hiZDCminus,fhi->hiZDCplus,(ptprobe-ptbarrel)/ptaverage,(float)(fskim->HBHENoiseFilterResultRun2Loose),(float)(fskim->pVertexFilterCutGplus),(float)(fhlt->HLT_AK4PFJet40_Eta5p1_v1),(float)(fhlt->HLT_AK4PFJet40_Eta5p1_v1_Prescl*fhlt->L1_SingleJet28_BptxAND_Prescl*fhlt->L1_SingleJet28_BptxAND_Prescl),(float)(fhlt->HLT_AK4PFJet60_Eta5p1_v1*fhlt->L1_SingleJet40_BptxAND_Prescl),(float)(fhlt->HLT_AK4PFJet60_Eta5p1_v1_Prescl*fhlt->L1_SingleJet40_BptxAND_Prescl),(float)(fhlt->HLT_AK4PFJet80_Eta5p1_v1),(float)(fhlt->HLT_AK4PFJet100_Eta5p1_v1),fhi->vz,fhi->vx,fhi->vy,(float)hasPForwardjet,(float)hasMForwardjet,randnum};
  float dijetEntry_corr[] = {pt1w,eta1w,phi1w,pt2w,eta2w,phi2w,pt3w,eta3w,phi3w,dphiw,trkmax1w,trkmax2w,ptaveragew,ptbarrelw,ptprobew,etabarrelw,etaprobew,pt3w/ptaveragew,fhi->hiHFplusEta4 + fhi->hiHFminusEta4,fhi->hiHFplusEta4,fhi->hiHFminusEta4,fhi->hiZDCminus,fhi->hiZDCplus,(ptprobew-ptbarrelw)/ptaveragew,(float)(fskim->HBHENoiseFilterResultRun2Loose),(float)(fskim->pVertexFilterCutGplus),(float)(fhlt->HLT_AK4PFJet40_Eta5p1_v1),(float)(fhlt->HLT_AK4PFJet40_Eta5p1_v1_Prescl*fhlt->L1_SingleJet28_BptxAND_Prescl),(float)(fhlt->HLT_AK4PFJet60_Eta5p1_v1),(float)(fhlt->HLT_AK4PFJet60_Eta5p1_v1_Prescl*fhlt->L1_SingleJet40_BptxAND_Prescl),(float)(fhlt->HLT_AK4PFJet80_Eta5p1_v1),(float)(fhlt->HLT_AK4PFJet100_Eta5p1_v1),fhi->vz,fhi->vx,fhi->vy,(float)hasPForwardjet,(float)hasMForwardjet,randnumw,pt1wup,pt1wdown,pt2wup,pt2wdown};
  ntdijet->Fill(dijetEntry);
  ntdijet_corr->Fill(dijetEntry_corr);
  ntdebug->Fill(pt1,eta1,trkmax1,trkSum1,chargedN1,rawpt1,eSum1);

}

ntdijet->Write();
ntjet->Write();
ntdijet_corr->Write();
ntdebug->Write();
fnt->Close();
}


int main(int argc, char *argv[])
{ 
  relativeResponse(argv[1],argv[2],argv[3]);
  return 0;
}