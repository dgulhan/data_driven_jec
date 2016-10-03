#include "L2L3ResidualWFits.h"
#include "MCTruthResidual.h"
#include "TH2D.h"
#include "TF1.h"

void test(){
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2(true);
  //smearing example
  TH2D * hist = new TH2D("hist","",20,50,400,20,0.8,1.2);
  L2ResidualJER* smear = new L2ResidualJER("pp5");
  cout << 1 << endl;
  TF1* fitResol = new TF1("fitResol","gaus(0)");
  fitResol->SetParameters(1.,1,0.2);
  double eta = 0;
  for(int i = 0; i < 100; i++){
   for(int ibin = 0; ibin < 20; ibin++){
    double pt = hist->GetXaxis()->GetBinCenter(ibin+1)*fitResol->GetRandom();
	double smearedPt = smear->getSmearedPt(pt, eta);
	hist->Fill(pt,smearedPt/hist->GetXaxis()->GetBinCenter(ibin+1));
   }
  }
  TFile *outf = new TFile("test.root","recreate");
  hist->Write();
  outf->Close();
  
  //pt correction example
  int radius = 3;
  double etacutcorr = 3;
  TString mode = "pp5";
  L2ResidualJES * L2JES = new L2ResidualJES(radius,((int)etacutcorr),mode);
  L2ResidualJER * L2JER = new L2ResidualJER(mode);
  L3ResidualJES * L3JES = new L3ResidualJES(mode);
  MCTruthResidual * MCTruth = new MCTruthResidual(mode); 

  
  double rawpt = 80.;
  double jtpt = 100.;
  double jteta = 0.;
  double correctedPt = jtpt;
  bool doMC = false;
  if(mode == "pPb5" || mode == "Pbp5") correctedPt = MCTruth->getJEC_1st(correctedPt,rawpt,jteta); 
  correctedPt = MCTruth->getResidualCorr(correctedPt,jteta);
  cout << "after MC correction jtpt = " << correctedPt << endl;
  if(!doMC){
   correctedPt = L2JES->getCorrectedPt(correctedPt, jteta);
   correctedPt = L3JES->getCorrectedPt(correctedPt);
   cout << "after data driven correction jtpt = " << correctedPt << endl;
  }else{   
   correctedPt = L2JER->getSmearedPt(correctedPt, jteta);
   cout << "after additional smearing for MC jtpt = " << correctedPt << endl;
  }
  
}