function [mag_out, phase_out] = c_spectral_mix(maga, magb, phasea, phaseb, gaina, gainb)
% function [mag_out, phase_out] = c_spectral_mix(maga, magb, phasea, phaseb, gaina, gainb)
%
% It computes the sum of two complex spectrograms with individual gains per channel.

c_out = (gaina* maga .* exp(1i*phasea))  + (gainb * magb .* exp(1i*phaseb));
mag_out = abs(c_out);
phase_out = angle(c_out);
end
