#!/bin/bash
# steering script to run the analysis code in correct order

echo "reading trees"
root -l -q src/readTaskOutput.C

echo "subtracting feed-down background"
root -l -q src/subtractLikeSigned.C

echo "reweighting"
root -l -q src/calculateWeights.C
root -l -q src/reweightReco.C
root -l -q src/findOptimalR.C

echo "calculating and applying AxE corrections"
root -l -q src/reweightAxE.C
root -l -q src/correctAxE.C

echo "getting omega phase parameter"
root -l -q src/fixOmega.C

echo "fitting"
root -l -q src/fitSpectra.C

echo "all done!"
