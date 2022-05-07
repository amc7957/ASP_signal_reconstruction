function SER = SER(original_signal, reconstructed_signal)
% RTISI_TSM(signal_in,fs) implements time scale modification using
% the RTISI algorithm presented in 'Real-Time Signal Estimation From 
%Modified Short-Time Fourier Transform Magnitude Spectra' by Xinglei Zhu, 
%Gerald Beauregard, and Lonce L. Wyse

os = original_signal;
rs = reconstructed_signal;
OS = fft(os,2048);
RS = fft(rs,2048);
SER = 10*log10((sum((1/2*pi)*abs(OS)).^2)/(sum((1/2*pi)*(abs(OS)-abs(RS)).^2)));
