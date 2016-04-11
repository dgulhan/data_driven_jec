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
  
  L2L3Residual(int radius=3, int etacut=3)
  {
   reset();
   this->radius = radius;
   this->etacut = etacut;
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
	
    /*
    eta_min[0] = -3;
    eta_max[0] = eta_min[1] = -2.500;
    eta_max[1] = eta_min[2] = -2.172;
    eta_max[2] = eta_min[3] = -1.740;
    eta_max[3] = eta_min[4] = -1.392;
    eta_max[4] = eta_min[5] = -1.044;
    eta_max[5] = eta_min[6] = -0.696;
    eta_max[6] = eta_min[7] = -0.348;
    eta_max[7] = eta_min[8] = 0.000;
    eta_max[8] = eta_min[9] = 0.348;
    eta_max[9] = eta_min[10] = 0.696;
    eta_max[10] = eta_min[11] = 1.044;
    eta_max[11] = eta_min[12] = 1.392;
    eta_max[12] = eta_min[13] = 1.740;
    eta_max[13] = eta_min[14] = 2.172;
    eta_max[14] = eta_min[15] = 2.500;
    eta_max[15] = 3;
	*/
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
    // eta_min[0] = -4;
    // eta_max[0] = eta_min[1] = -3.664;
    // eta_max[1] = eta_min[2] = -3.314;
    // eta_max[2] = eta_min[3] = -2.964;
    // eta_max[3] = eta_min[4] = -2.650;
    // eta_max[4] = eta_min[5] = -2.322;
    // eta_max[5] = eta_min[6] = -2.043;
    // eta_max[6] = eta_min[7] = -1.830;
    // eta_max[7] = eta_min[8] = -1.653;
    // eta_max[8] = eta_min[9] = -1.479;
    // eta_max[9] = eta_min[10] = -1.305;
    // eta_max[10] = eta_min[11] = -1.131;
    // eta_max[11] = eta_min[12] = -0.957;
    // eta_max[12] = eta_min[13] = -0.783;
    // eta_max[13] = eta_min[14] = -0.609;
    // eta_max[14] = eta_min[15] = -0.435;
    // eta_max[15] = eta_min[16] = -0.261;
    // eta_max[16] = eta_min[17] = -0.087;
    // eta_max[17] = eta_min[18] = 0.087;
    // eta_max[18] = eta_min[19] = 0.261;
    // eta_max[19] = eta_min[20] = 0.435;
    // eta_max[20] = eta_min[21] = 0.609;
    // eta_max[21] = eta_min[22] = 0.783;
    // eta_max[22] = eta_min[23] = 0.957;
    // eta_max[23] = eta_min[24] = 1.131;
    // eta_max[24] = eta_min[25] = 1.305;
    // eta_max[25] = eta_min[26] = 1.479;
    // eta_max[26] = eta_min[27] = 1.653;
    // eta_max[27] = eta_min[28] = 1.830;
    // eta_max[28] = eta_min[29] = 2.043;
    // eta_max[29] = eta_min[30] = 2.322;
    // eta_max[30] = eta_min[31] = 2.650;
    // eta_max[31] = eta_min[32] = 2.964;
    // eta_max[32] = eta_min[33] = 3.314;
    // eta_max[33] = eta_min[34] =  3.664;
    // eta_max[34] = 4;
   }
   correction_file = new TFile(Form("L2L3VsPtEtaBinned_alphacut_high2_%s_etacut%d.root",algo_corr.Data(),etacut));
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