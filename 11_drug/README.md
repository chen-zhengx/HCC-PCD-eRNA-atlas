# 11_drug: drug sensitivity prediction and candidate drugs

This module predicts the response (IC50) of 198 GDSC2 drugs in the 70 paired tumor and adjacent non-tumor samples of
GSE144269 with oncoPredict, tests for each drug whether the predicted IC50 is lower in the tumor than in the paired
adjacent tissue, and reproduces the 19 candidate drugs of Table S11 together with the full results for all 198 drugs.
