#include <cmath>
#include <iostream>
#include "TF1.h"
#include "TFile.h"

using namespace std;
class smearDataMCJER
{
  private:
    TF1 *fitSF[2];
    TF1 *fitGaus;
	
    const double ptLow = 50;
    const double ptHigh = 400;
	const double etaBound = 1.3;
	const double etaMax = 3.;
	int neta = 2;
	TString mode; //pp or pPb
  public:
  smearDataMCJER(TString mode){
    
	double parsPPb[2][2] = {{5.61139e+00,9.40866e-01},
                       {1.19898e+00,6.28925e-01}};
	  
	double parsPP[2][2] = {{7.57848e-01,6.23244e-01},
	                  {6.08142e+00,1.05733e+00}};
	  
	for(int ieta = 0; ieta < neta; ieta++){
     fitSF[ieta] = new TF1(Form("fitSF_eta%d",ieta),"[0]/pow(x,[1])",50,400);
    }
	
	this->mode = mode;
	for(int ieta = 0; ieta < neta; ieta++){
	  for(int ipar = 0; ipar < 2; ipar++){
  	    if( mode == "pPb" ){
		  fitSF[ieta]->SetParameter(ipar, parsPPb[ieta][ipar]);
		}else{
		  fitSF[ieta]->SetParameter(ipar, parsPP[ieta][ipar]);
		}
	  }
	}
    fitGaus = new TF1("fitGaus","gaus(0)",-20,20);
    fitGaus->SetParameters(1,0,1);
  }
 
  int getIEta(double eta){
   if( fabs(eta) < 1.3 ) return 0;
   else return 1;
  }
 
 double getSmearedPt(double pt, double eta){
   if( pt < ptLow || pt > ptHigh || fabs(eta) > etaMax ) return pt;
   cout << "passed" << endl;
   cout << eta << " " << getIEta(eta) << endl;
   cout << fitGaus->GetRandom()<<endl;
   return pt*(1+(fitSF[getIEta(eta)]->Eval(pt))*fitGaus->GetRandom());
 }
};