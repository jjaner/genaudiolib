function pitch_cont = get_continuous_pitch(pitch)
% function pitch_cont = get_continuous_pitch(pitch)
%
% Interpolates pitch values for unvoiced frames (pitch==0).
%
% Input: 
% - pitch: array with one pitch per frame in Hz. Unvoiced have pitch(i) =
% 0.
%
% Output: 
% - pitch_cont: output array with one pitch per frame in Hz.

pitch_cont = pitch;
unvcB = find((pitch(2:end)==0) .* (diff(pitch) <0) );
unvcE = find((pitch(1:end-1)==0) .* (diff(pitch) >0) );
nFrames = size(pitch,1);
pitch(isnan(pitch)) = 0;
% plot(pitch);
% hold on
% stem(unvcB,ones(size(unvcB)),'rs');
% stem(unvcE,ones(size(unvcE)),'gx');

%for i=1:size(unvcB)
if pitch(1) == 0
    unvcB = [1; unvcB];
end
if pitch(end) == 0
    unvcE = [unvcE; nFrames];
end

itype = 'linear';
for i= 1:size(unvcB)
   
    leftP = pitch(unvcB(i));
    if  unvcB(i)== 1 && pitch(unvcB(i)) == 0,
            leftP = pitch(min(nFrames,unvcE(i)+1)); % first frame is unvoiced 
    end
    rightP = leftP; 
    if unvcE(i) == nFrames,
        rightP = pitch(unvcB(i)); 
    else
        rightP = pitch(unvcE(i)+1); 
    end    
    pint  = interp1([unvcB(i), unvcE(i)+1],[leftP,rightP], [unvcB(i)+1:unvcE(i)],itype);
    pitch_cont(unvcB(i)+1:unvcE(i)) = pint;
end


end
          