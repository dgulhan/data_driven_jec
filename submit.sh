#!/bin/sh


version=1
#dataset=/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories//2015pp_HighPtJet80
dataset=/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_HighPtLowerJets/
#dataset=/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_MinBias6
# dataset=/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_MinBias_2
destination=/mnt/hadoop/cms/store/user/dgulhan/ppDijet/20160408/lowerpt_ak3PF_eta3_negshift

NAME="relativeResponse.C"
# NAME="relativeResponse_MC.C"
g++ $NAME $(root-config --cflags --libs) -O2 -o "${NAME/%.C/}.exe"

subdir=`date +%Y%m%d`
logdir=/net/hisrv0001/home/dgulhan/logsDijet/${subdir}_v3
mkdir $logdir

mkdir -p $destination
# cp Corrections/Casym_pp_double_hcalbins_algo_ak4PF_pt100_140_jet80_alphahigh_20_phicut250.root .
cp Corrections/L2L3VsPtEtaBinned_*.root .
cp Corrections/someCorr/* .
# tar -cvzf corr.tar.gz trkCorrections

for file in ${dataset}/*
  do
      cropped=`echo $file | sed "s/\// /g" | awk '{ print $9 }'`
	  echo $cropped
	  outfile=hltTree_${cropped}
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
Executable     = runTreeMaker.sh

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
transfer_input_files    = relativeResponse.exe,L2L3VsPtEtaBinned_alphacut_high2_ak3PF_etacut4.root,Casym_pp_hcalbins4_algo_ak3PF_pt75_100_jet60_corrv5_eta4_alphahigh_20_phicut250_etacut4.root,Casym_pp_hcalbins4_algo_ak3PF_pt55_75_lowerpt_corrv5_eta4_alphahigh_20_phicut250_etacut4.root,Casym_pp_hcalbins4_algo_ak3PF_pt100_400_jet80_corrv5_eta4_alphahigh_20_phicut250_etacut4.root

Queue
EOF

# submit the job
condor_submit subfile

done

# rm Casym_pp_double_hcalbins_algo_ak4PF_pt100_140_jet80_alphahigh_20_phicut250.root


