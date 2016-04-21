#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

using namespace std;
class L2L3Residual
{
 private:
  int neta;
  int etacut;
  double lower_pt_cut;
  double higher_pt_cut;
  int radius;
  bool dopPb;
  TFile *correction_file;
  TF1 * fits[100];
  TString algo_corr;
  double eta_min[100];
  double eta_max[100];
  public:
  
   void reset()
   { 
    for(int ieta=0;ieta<100;ieta++){
     fits[ieta] = NULL;
    }

    correction_file = NULL;
    
   }
  
  L2L3Residual(int radius=3, int etacut=3, bool dopPb=false)
  {
   reset();
   this->radius = radius;
   this->etacut = etacut;
   this->dopPb = dopPb;
   algo_corr = Form("ak%dPF",radius);   
   
   lower_pt_cut = 30;
   higher_pt_cut = 400;
   
   if(etacut==3){
    // neta=16;
    neta=58;
    double etabins_hcalbins[]= {-3, -2.853,
                     -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830,
                     -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
                     -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609,
                     -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,  0.000,
                      0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,
                      0.696,  0.783,  0.879,  0.957,  1.044,  1.131,  1.218,
                      1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,
                      1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  2.853,
                      3};
					  
    for(int ieta = 0; ieta < neta; ieta++){
	 eta_min[ieta] = etabins_hcalbins[ieta];
	 eta_max[ieta] = etabins_hcalbins[ieta+1];
	}
	
   }
   if(etacut==4){
	neta=64;
	double etabins_hcalbins4[]= {-4,      -3.664,  -3.314,  -2.964, -2.853,
                     -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830,
                     -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
                     -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609,
                     -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,  0.000,
                      0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,
                      0.696,  0.783,  0.879,  0.957,  1.044,  1.131,  1.218,
                      1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,
                      1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  2.853,
                      2.964,  3.314,  3.664, 4};
    for(int ieta = 0; ieta < neta; ieta++){
	 eta_min[ieta] = etabins_hcalbins4[ieta];
	 eta_max[ieta] = etabins_hcalbins4[ieta+1];
	}
   }
   correction_file = new TFile(Form("L2L3VsPtEtaBinned_alphacut_high2_%s_etacut%d_dopPb%d.root",algo_corr.Data(),etacut,dopPb));
   for(int ieta = 0;ieta<neta;ieta++){
    fits[ieta] = (TF1*)correction_file->Get(Form("fit%d",ieta));
   } 
  }
  
  double get_corrected_pt(double jetpt, double jeteta)
  {
   double correction = 1;
   if( abs(jeteta)> ((double)etacut)) return correction*jetpt;
   int etaindex = 0;
   for(int ieta = 0; ieta < neta; ieta++){
    if(eta_min[ieta] > jeteta ) continue;
	else etaindex = ieta;
   }

   if(jetpt < lower_pt_cut) return fits[etaindex]->Eval(lower_pt_cut)*jetpt;
   if(jetpt > higher_pt_cut) return fits[etaindex]->Eval(higher_pt_cut)*jetpt;

   return fits[etaindex]->Eval(jetpt)*jetpt;
  
  }
};