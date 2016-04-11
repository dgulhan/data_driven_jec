#!/bin/sh


version=1
#dataset=/mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet15_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_125350/0000/
# dataset=/mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet30_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133336/0000/
# dataset=/mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet50_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133347/0000/
#dataset=/mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet80_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133359/0000/
#  dataset=/mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet120_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133414/0000/
dataset=/mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet170_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133425/0000/
# dataset=/mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet220_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133441/0000/
# dataset=/mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet280_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133453/0000/
# dataset=/mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet370_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133505/0000/
destination=/mnt/hadoop/cms/store/user/dgulhan/ppDijet/20160401/pthat170_ak3PF_eta3_shift

# NAME="relativeResponse.C"
NAME="relativeResponse_MC.C"
g++ $NAME $(root-config --cflags --libs) -O2 -o "${NAME/%.C/}.exe"

subdir=`date +%Y%m%d`
logdir=/net/hisrv0001/home/dgulhan/logsDijet/${subdir}_v3
mkdir $logdir

mkdir -p $destination

# tar -cvzf corr.tar.gz trkCorrections

for file in ${dataset}/*
  do
      cropped=`echo $file | sed "s/\// /g" | awk '{ print $11 }'`
	  echo $cropped
	  outfile=hltTree_HighPtLowerJets_${cropped}
      Error=`echo $outfile | sed "s/root/err/g"`
	  Output=`echo $outfile | sed "s/root/out/g"`
	  Log=`echo $outfile | sed "s/root/log/g"`        
	  	
	      echo "Output is : $outfile"
	      echo "Error is : $Error"
	      echo "LFN is : $lfn"
	      echo "----------------------"
	      
	      cat > subfile <<EOF

Universe       = vanilla

# files will be copied back to this dir
Initialdir     = .

# run my script
Executable     = runTreeMakerMC.sh

+AccountingGroup = "group_cmshi.yetkin"
#+IsMadgraph = 1

Arguments      = $dataset $cropped $destination $outfile \$(Process)
# input files. in this case, there are none.
Input          = /dev/null

# log files
Error          = $logdir/$Error
Output         = $logdir/$Output
Log            = $logdir/$Log

# get the environment (path, etc.)
Getenv         = True

# prefer to run on fast computers
Rank           = kflops


# should write all output & logs to a local directory
# and then transfer it back to Initialdir on completion
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
# specify any extra input files (for example, an orcarc file)
transfer_input_files    = relativeResponse_MC.exe

Queue
EOF

# submit the job
condor_submit subfile

done



