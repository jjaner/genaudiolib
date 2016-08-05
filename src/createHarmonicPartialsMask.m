
function mask = createHarmonicPartialsMask(pitch, harmonics_list, param)

specsize = (param.fftsize/2)+1;
freq2bin = 2*specsize/param.fs;

nFr = size(pitch,1);
mask = zeros(nFr, specsize);
width = 75; % percentage of distance between partials

for i=1:nFr,
    listIdx = find(harmonics_list(:,1)<i);
    listIdx = listIdx (end);
    bypassHarms = [harmonics_list(listIdx,2):harmonics_list(listIdx,3)];
    if pitch(i)>35 % 60Hz
        % set to ones at positions of the by-passed harmonics      
        harmIdx = floor([pitch(i): pitch(i):param.fs/2]*freq2bin); 
        bypassHarms = min(bypassHarms, size(harmIdx,2));
        partialWidth = floor((width/100) * pitch(i) *  freq2bin);
        mask(i,min(specsize,harmIdx(bypassHarms))) = 1;  
        h = [ones(1,partialWidth)];
        % mask will set to ones partialWidth around the harnmonic
        mask(i,:) = conv(mask(i,:),h,'same');
    end
end
    % smooth output mask over time for softer transitions
%     widthTemporal = 5;
%     widthSpectral = 20;
%     sig = [widthTemporal widthSpectral]; % sigma(1): temporal in frames, sigma(2): spectral in bins
%     N = max(sig)*5;
%     ker = gaussian2d(N,sig);
%     ker = ker(:,round(size(ker,1)/2)- sig(1)*2:round(size(ker,1)/2) +sig(1)*2);
%     ker = ker/max(max(ker));
%     mask = conv2(mask,ker,'same'); % smooth    
    h = [0.1 0.25 1 0.25 0.1];
    
    %original function
    %mask = filter2(h/sum(h),mask')'; % smooth
    % alternative reduce memory requirement
    specSize = size(mask,2);
    specStep  = floor(specSize/128);    
    for j=1:specStep:specSize-specStep,        
        mask(:,j:j+specStep) = filter2(h/sum(h),mask(:,j:j+specStep)')';
    end
end


% -------------------------------------------
% ------------------------------------
function f=gaussian2d(N,sigma)
  % N is grid size, sigma speaks for itself
  
 [x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
 f=exp(-x.^2/(2*sigma(1)^2)-y.^2/(2*sigma(2)^2));
 f=f./sum(f(:));
 end
