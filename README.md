# Kayseretal_EJN24_ccj

Data and analysis code  for the publication 
"Perceived multisensory common cause relations shape the ventriloquism effect but only marginally the trial-wise aftereffect" 
by Christoph Kayser, Herbert Heuer
Department of Cognitive Neuroscience, Universität Bielefeld, Bielefeld, Germany Leibniz Research Centre for Working Environment and Human Factors, Dortmund, Germany Contact: Christoph.kayser@uni-bielefeld.de.

This study featured two experiments consisting of sequences of audio-visual and auditory trials in which we quantified the ventriloquist effect and the trial-wise aftereffect. The main feature is the inclusion of common-cause judgements in the AV trials, contingent on which we analyze the data. The main hypothesis is whetehr each effect depends on the common-cause judgement. 

The folder contains the preprocessed data files for each experiment:
CK12_Exp1_Alldata.mat and CK12_Exp2_Alldata.mat

These contain a cell array with the data for each participant AllData{participant}
this features for each pair of AV-A trials the following information

AllData{participant}(:, [VEbias, VAEbias, Vpos, Apos, DVA, Rel, APosA, respAV, respA, CCJ, Trial  sub]).

Here VEbias is the ventriloquist bias obtained in the AV trial, VAEbias the aftereffect in the A trial, 
Vpos and Apos the stimulus positions in the AV trial, DVA their difference (Vpos-Apos). Rel is the reliability 
of the acoustic stimulus (1 for high, 2 for low). AposA the sound position in the A trial. respAV (respA) the continous response 
in the AV (A) trial. CCj the binary confidence judgement for the AV trial (2 yes, 1 no) . Trial and Sub the trial and participant index.


Analyze_CK12_Analysis_IndivSlopes.m performs the main analysis for each experiment. 
Create the figure, computes single participant bias slopes and runs an lme testing the effect of reliability and ccj on this. 
Some parts require additional matlab packages (see comments in code).
