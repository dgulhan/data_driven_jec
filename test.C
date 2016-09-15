#include "smearDataMCJER.h"
#include "TH2D.h";

void test(){
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2(true);
  
  TH2D * hist = new TH2D("hist","",20,50,400,20,0.8,1.2);
  cout<<0<<endl;
  smearDataMCJER* smear = new smearDataMCJER("pp");
  cout << 1 << endl;
  TF1* fitResol = new TF1("fitResol","gaus(0)");
  fitResol->SetParameters(1.,1,0.2);
  for(int i = 0; i < 100; i++){
  for(int ibin = 0; ibin < 20; ibin++){
    cout<<2<<endl;
    double x = hist->GetXaxis()->GetBinCenter(ibin+1)*fitResol->GetRandom();
	double smearedPt = smear->getSmearedPt(x, 0);
	hist->Fill(x,smearedPt/x);
	cout<< 3<< endl;
  }
  }
  TFile *outf = new TFile("test.root","recreate");
  hist->Write();
  outf->Close();
  
}