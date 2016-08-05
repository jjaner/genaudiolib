function [fo, mo] = midiFreqConversion(mi,fi)
% function [fo mo] = midiFreqConversion(mi,fi)
% MIDI to Frequency bijective conversions 
% 
% mi: MIDI note input
% fi: frequency input in Herz
%
% mo: MIDI note output
% fo: frequency output in Herz

if size(mi)>0
    fo = (440/32)*2.^((mi-9)/12);
    mo = mi;
end
if size(fi)>0
    mo = 9 + 12 * log(fi / (440/32))/(log(2)) ;
    fo = fi;
end
