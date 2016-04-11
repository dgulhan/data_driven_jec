# Jet energy correction files

Data driven JEC calculation relative correction with dijet balance and absolute correction with gamma-jet balance
Steps:
1. Generate ntuples for data, each jet trigger separately, and MC, each pthat separately. Only a single pthat sample and trigger is used in a given pt bin in derivation. 
For data:
```
 g++ relativeResponse.C $(root-config --cflags --libs) -O2 -o relativeResponse.exe
 ./relativeResponse.exe <directoryname> <filename> <outputfilename>
```
For MC:
```
 g++ relativeResponse_MC.C $(root-config --cflags --libs) -O2 -o relativeResponse_MC.exe
  ./relativeResponse_MC.exe <directoryname> <filename> <outputfilename>
```
I was using bash scripts with condor to do this.

While running the above macros make sure the following parameters match in the relativeResponse.C and relativeResponse_MC.C:
 float etacut=3;
 TString algo="ak3PF";
 
2. Edit the file calculateCasymVsPtEtaBinned.C to have the merged output of the ntupler.
3. Run calculateCasymVsPtEtaBinned.C. In Corrections folder root files will appear
4. Check the fits and balance distributions by plotting them using plotCorrections/plotBandMean.C
5. Produce ntuples from photon forests to derive absolute corrections by running absoluteResponse.C both for data and MC doing the same as for relativeResponse.C
6. Run extraAlpha.C 
7. Run calculateAbsRes.C
