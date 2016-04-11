for ifile in 0
 do
 output=out_${file[$ifile]}

  echo | awk -v dataset=$1 -v infile=$2 -v outfile=$4 '{print "./relativeResponse_MC.exe \""dataset"\" \""infile"\" \""outfile"\""}' | bash

  mv $4 $3
 done
