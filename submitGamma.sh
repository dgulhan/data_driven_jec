#!/bin/sh


version=1
# dataset=/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories//2015pp_HighPtJet80
# dataset=/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_HighPtLowerJets/
# dataset=/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_MinBias6
#dataset=/mnt/hadoop/cms/store/user/richard/2015-Data-photonSkims/HighPtPhoton30AndZ/pp-photonSkim-wFilter-v0/160106_185221/0000/ #photon30
# dataset=/mnt/hadoop/cms/store/user/richard/2015-Data-promptRECO-photonSkims/HighPtPhoton30AndZ/pp-photonHLTFilter-v0/160125_194718/0000/ #photon40
# dataset=/mnt/hadoop/cms/store/user/rbi/Pythia8_Photon30_pp502_TuneCUETP8M1/Pythia8_Photon30_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_forest_v1/160315_201109/0000/ #pthat30
#dataset=/mnt/hadoop/cms/store/user/rbi/Pythia8_Photon15_pp502_TuneCUETP8M1/Pythia8_Photon15_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_forest_v1/160315_200126/0000/ #pthat15
#dataset=/mnt/hadoop/cms/store/user/rbi/Pythia8_Photon50_pp502_TuneCUETP8M1/Pythia8_Photon50_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_forest_v1/160315_201159/0000/
dataset=/mnt/hadoop/cms/store/user/rbi/Pythia8_Photon120_pp502_TuneCUETP8M1/Pythia8_Photon120_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_forest_v1/160315_201228/0000/
# dataset=/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_MinBias_2
# destination=/mnt/hadoop/cms/store/user/dgulhan/ppDijet/20160330/photon30_ak3PF
# destination=/mnt/hadoop/cms/store/user/dgulhan/ppDijet/20160330/photon40_ak3PF
destination=/mnt/hadoop/cms/store/user/dgulhan/ppDijet/20160330/pythia_pho120_ak3PF
# destination=/mnt/hadoop/cms/store/user/dgulhan/ppDijet/20160330/pythia_pho30_ak3PF

NAME="absoluteResponse.C"
# NAME="relativeResponse_MC.C"
g++ $NAME $(root-config --cflags --libs) -O2 -o "${NAME/%.C/}.exe"

subdir=`date +%Y%m%d`
logdir=/net/hisrv0001/home/dgulhan/logsDijet/${subdir}_v3
mkdir $logdir

mkdir -p $destination
# cp Corrections/Casym_pp_double_hcalbins_algo_ak4PF_pt100_140_jet80_alphahigh_20_phicut250.root .
cp Corrections/L2L3VsPtEtaBinned_ak3PF.root .
# tar -cvzf corr.tar.gz trkCorrections

for file in ${dataset}/*
  do
      cropped=`echo $file | sed "s/\// /g" | awk '{ print $11 }'`
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
Executable     = runTreeMakerGamma.sh

+AccountingGroup = "group_cmshi.dgulhan"
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
transfer_input_files    = absoluteResponse.exe,L2L3VsPtEtaBinned_ak3PF.root

Queue
EOF

# submit the job
condor_submit subfile

done

# rm Casym_pp_double_hcalbins_algo_ak4PF_pt100_140_jet80_alphahigh_20_phicut250.root


