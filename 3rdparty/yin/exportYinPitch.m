function dataout = exportYinPitch(wavpath, hopsize)

P.hop = hopsize; %256; %512; % s - interval between estimates (default: 32/SR)
P.maxf0 = 1800;

R = yin(wavpath,P);
outfile = strrep(wavpath,'wav','pitch');
outfile = strrep(outfile,'WAV','pitch');
outfileTools = strrep(outfile,'pitch','tools.f0');

pitchout  = [440*2.^(R.f0) .* (((R.ap0).^2) < 0.35)];
hopsizeRatio = 1; %6 % defulat yin hopsize is 32/SR. Ours is 512/44100

pitchout512 = pitchout(1:hopsizeRatio:end);
%time512 = [0.0639:512/44100:(size(pitchout512,2)-1)*512/44100];
time512 = [2*P.hop/44100:P.hop/44100:(size(pitchout512,2)-1)*P.hop/44100];
dyn = 127*R.pwr(1:hopsizeRatio:end)'/max(R.pwr);
dyn = dyn(1:size(time512,2));
pitchout512 = pitchout512(1:size(time512,2));
% remove NaN
dyn = dyn .* (dyn == dyn);
pitchout512(isnan(pitchout512)) = 0;
dyn(isnan(dyn)) = 0;

% first frames offset
initframes = [0 , 0 ,0 ; P.hop/44100,0,0];
% write output file
dataout = [initframes ; time512' dyn pitchout512'];
csvwrite(outfile,dataout);
save(outfileTools,'dataout','-ascii');