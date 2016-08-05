function [psAfact, psBfact, synTimeIdxA, synTimeIdxB, morphAB] = read_score_edwin(score_file, param)
% function [psAfact, psBfact, synTimeIdxA, synTimeIdxB, morphAB] = read_score_edwin(score_file, param)
%
% Score data is passed as: ot, i1t, i1p, i2t, i2p, mph
%  - Times are in seconds
%  - Pitch shift are in semitones (+/-)
%  - Morphing values are [0..1]
% Outputs:
%  - PSAfact, psBfact: pitch shifting factor arrays (one value per output frame)
%  specified in linear factors (e.g. psAfact=2 means +1 octave transposition)
%  - synTimeIdxA: index of the frame index in the input sound (one value per output frame). It
%  oversamples (time expansion) or sub-samples (tome compressions) the
%  frame index array.


data = load(score_file);
nPoints = size(data,1);

% TS file A
% variable TS - synthesis time index  (tine in time out in seconds)
% times_inout  =  [0 0 ;0.5 2.5 ;1.95 5.4; 2.6 6.05];
times_inout  = [data(:,2) data(:,1)]; % 
morph_in = data(:,6);
ps_in = data(:,3);
nFrA = data(nPoints,2);

synTimeIdx = []; % start inddexes at 1
morphAB = [];
psA = [];
for i=1:size(times_inout,1)-1,
    bfri = (times_inout(i,1) * param.fs / param.hopsize);
    bfr = (times_inout(i,2) * param.fs / param.hopsize);
    efr = (times_inout(i+1,2) * param.fs / param.hopsize);
     vect_new = [0:(efr-bfr)-1];
    ts_tmp = (times_inout(i+1,2) - times_inout(i,2)) / (times_inout(i+1,1) - times_inout(i,1));
    synTimeIdx = [synTimeIdx  , 1 + bfri + [(0:(efr-1)-bfr)/ts_tmp]];
    
    % pitch shift
    bps = ps_in(i);
    eps = ps_in(i+1);  
    psA = [psA  , bps + ((eps-bps)/(efr-bfr)) * vect_new];
    
    % morph
    bmph = morph_in(i);
    emph = morph_in(i+1);   
    morphAB = [morphAB  , bmph +  ((emph-bmph)/(efr-bfr))*vect_new];

end
synTimeIdxA = synTimeIdx';
morphAB = morphAB';

% TS file B
times_inout  = [data(:,4) data(:,1)]; %  [0 0 ;0.5 2.5 ;1.95 5.4; 2.6 6.05];
nFrB = data(nPoints,4);
ps_in = data(:,5);

synTimeIdx = []; % start inddexes at 1
psB = [];
for i=1:size(times_inout,1)-1,
    bfri = (times_inout(i,1) * param.fs / param.hopsize);
    bfr = (times_inout(i,2) * param.fs / param.hopsize);
    efr = (times_inout(i+1,2) * param.fs / param.hopsize);
     vect_new = [0:(efr-bfr)-1];
    ts_tmp = (times_inout(i+1,2) - times_inout(i,2)) / (times_inout(i+1,1) - times_inout(i,1));
    synTimeIdx = [synTimeIdx  , 1 + bfri + [(0:(efr-1)-bfr)/ts_tmp]];
    
    % pitch shift
    bps = ps_in(i);
    eps = ps_in(i+1);   
    psB = [psB  , bps + ((eps-bps)/(efr-bfr)) * vect_new];

end
synTimeIdxB = synTimeIdx';

% convert pitch shifting from semitones to linear pitch shifting factors
psAfact = 2.^(psA/12);
psBfact = 2.^(psB/12);

end