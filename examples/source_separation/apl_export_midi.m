function apl_export_midi(out_midi, outmidifile)
% function apl_export_midi(out_midi, outmidifile)
%
% Input:
%   out_midi: matrix with MIDI format as specified in the MIDILib
% Output:
%   outmidifile: output MIDI filename (.mid)
%
% Example function of MIDI file export from matlab data
%-------------------------------------------------------------
% beat info ticks per quarter note

% MIDI file parameters
param.midi_tpqn = 120;
param.midi_bpm = 120;


if size(out_midi)>0
  [out_folder ~] = fileparts(outputMIDIfile);
  if ~exist(out_folder),
      disp('Warning: MIDI file was not created because the target folder does not exist.')
else
  writemidi_java(out_midi, outmidifile, param.midi_tpqn, param.midi_bpm);
  end
else
  disp('Warning: MIDI file was not created because the note transcription was empty.')
end

end  % end function

