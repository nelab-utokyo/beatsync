# code

## shuf_param2.py

This code is a generator of .mat file from .txt file saved for the experimental parameters.
TC, and stimuli order in the experiment is converted.
shuf_param.txt -> shuf_param.mat


## data_extract_ver5.m

LFP data and MUA data is extracted from the recorded data.
ns2, nev2 -> exdata_ex2.mat etc.


## evoked_response.m

This code is to produce the auditory response MUA.
data2matrix3 function is used in this code.
exdata_ex2.mat -> evoke_res3.mat


# Procedure

1. Implement TonotopicMap_ALL_uta2016 to get Characteristic frequency of each electrode.
2. Manually make shuffle_parameters.txt and write down the electrode number with tune curve, and music order.
3. Implement shuf_param2.py and obtain shuf_param.mat.
4. Implement data_extract_ver5.m to trim the related data and save it in .mat file.
5. Implement evoked_response.m to calcurate the evoked MUA.
