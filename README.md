

DeepTauID
==============


git clone https://github.com/jkiesele/DeepTauID DeepTau


The main 'cmsRun' script is here: 
https://github.com/jkiesele/DeepTauID/blob/master/DeepTauID/production/DeepTauNTuples.py

Please ignore the 'test' directory.

There is a job submission script for (lxplus) condor here:
https://github.com/jkiesele/DeepTauID/blob/master/DeepTauID/scripts/jobSub.py

Coming with a check script:
https://github.com/jkiesele/DeepTauID/blob/master/DeepTauID/scripts/check.py

THe actual defininition of the modules filling the variables can be found in the ntuple_X files here:
https://github.com/jkiesele/DeepTauID/tree/master/DeepTauID/src


The file:
https://github.com/jkiesele/DeepTauID/blob/master/DeepTauID/plugins/DeepTauNTuplizer.cc
contains a few globals and the gen matching parts. For filling the ntuple, it just calls the modules mentioned above.
