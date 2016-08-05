 -------- MIESPV CODE ------------------------------------------------
 MIESPV - Monophonic audIo procESsing based on Phase Vocoder
 ---------------------------------------------------------------------
 This software implements a set of sound transformations based on phase vocoder technique. The provided sound transformations are adaptations of extiting techniques developed at the MTG-UPF, adapted to meet the artistic goals
 of the sound installation authored by Edwin van der Heide in the context of SONAR 2014 in Barcelona. 

 The MIESPV code consists of three components:
 - MiesPV.m: A Matlab script adapted to produce the sound transformations required for the Edwin van der Heide's installation. It loads a score file, and outputs a number of transformed audio files for all harmonic partials. This script is developed by Jordi Janer and  Antonio Sá. The code is propierty of MTG-UPF.
 - GenaudioLib library: a Matlab code library (GenAudioLib) for STFT analysis/synthesis and Phase Vocoding. This library is proprierty of the MTG-UPF and developed/maintained by Jordi Janer. Provided as pre-compiled code. 
 - YIN library: a third-party Matlab library (YIN) for monophonic pitch estimation. This library is NOT propierty of the MTG-UPF.

 
 Collaboration of Music Technology Group/Universitat Pompeu Fabra, Sonar+D
 2014 and Edwin van der Heide.

 SMC student: Antonio Sá
 MTG supervisor: Jordi Janer
 Barcelona May 2014.



USAGE:
>> MiesPV


To run the code, you should first modify the section USER PARAMETERS in the script MiesPV.m.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER PARAMETERS

% input / output files
inputFileName = '../edwin/erfan.wav';
outfile = '../edwin/out/erfan_pv_TS_2_PS_2.wav';

% Test input duration
testDuration = 3; % limits the duration of the input sound to the specified value in seconds. Default: -1 (original size is used)

% Transformation parameters
TS = 2; % in time-scaling factor = dur_out/ dur_in. Default: 1
PS = 2;  % in semitones. Default: 0
nHarmonics = 40; % number of output harmonic partials (one file per harmonic will be created)
useTimbrePreserve = false; % Use timbre (formant) perservation. Set to true for large pitch shift values. Default: false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



