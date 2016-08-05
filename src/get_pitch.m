function [pitch_out midi_mat] = get_pitch(pitchfilename, nframes, hopsize, srate)

pitch_out = [];
t_o = [];
midi_mat = [];

[pathstr,name,ext] = fileparts(pitchfilename);
isMIDI = strcmp(ext,'.mid') || strcmp(ext,'.MID');


if exist(pitchfilename)
    pitch = [];
    if (isMIDI) % TODO
        midi_mat = readmidi_java(pitchfilename);
        disp('TODO: get pitch informatiopn from MIDI')
        % TODO convert MIDI data to pitch per frame info...
        midi2Hz_path = 'D:\XSource\roma\monet\utils\python\midi2Hz\midi2hz.py';
        cmd = sprintf('python %s %s -a',midi2Hz_path, pitchfilename);
        %system(cmd);
        %pitchfilename = strrep(pitchfilename,'.mid','.csv');
    else
        pitch = load(pitchfilename,'-ascii');
        pitchnframes = size(pitch,1);
        pitchntracks = (size(pitch,2) -1)/2;
        
        t_o = pitch(:,1);
        if pitchntracks == 1,
            % one track pitch file format
            pitch = pitch(:,3);
        else
            % multi-track pitch file format
            pitch = pitch(:,2:2:end);
        end
        
        dur = t_o(end);
        nframes = floor(dur*srate / hopsize);
        pitch_out = zeros(size(pitch));
        pitch_out(~isnan(pitch)>0) = pitch (~isnan(pitch)>0);
        
        % resample to match the target hopsize
        t_t = [0:hopsize/srate:nframes*hopsize/srate]'; % target times
        pinter_near = interp1(t_o,pitch_out,t_t,'nearest');
        pinter_line = interp1(t_o,pitch_out,t_t,'linear');
        pitch_out = (pinter_near > 0) .* pinter_line;
        pitch_out = pitch_out(1:nframes,:);
    end % end isMIDI
    
else    
    sprintf('Computing pitch file %s with Yin.',pitchfilename);
    wavpath = strrep(pitchfilename,'.pitch','.wav');
    pitch_out = exportYinPitch(wavpath, hopsize);    
    pitch_out = pitch_out(:,3);
end
end
