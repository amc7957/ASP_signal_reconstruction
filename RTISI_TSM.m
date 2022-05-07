function RTISI_TSM(signal_in, fs, alpha)
% RTISI_TSM(signal_in,fs) implements time scale modification using
% the RTISI algorithm presented in 'Real-Time Signal Estimation From 
%Modified Short-Time Fourier Transform Magnitude Spectra' by Xinglei Zhu, 
%Gerald Beauregard, and Lonce L. Wyse. alpha = 2 makes the signal speed 
%up and alpha = 0.5 makes the signal slow down

[row,col] = size(signal_in);
%L window length
L = 256; 
%Ss step size
Ss = L/4; 
a = 0.54;
b = -0.46;
n = 1:L;
w = (2*sqrt(Ss))/sqrt((4*a^2+2*b^2)*L)*(a+b*cos(2*pi*n/L)).';
iterations = 5;
Sa = Ss/alpha; 
par_frame = zeros(row/alpha,1);

for n = 0:((row-L)/Ss)
    
    %create partial frame
    if L+(n*Sa)<=row/alpha
        temp_frame = par_frame(n*Sa+1:L+(n*Sa));
    %make last frame in signal the right frame size    
    else
        temp_frame = [par_frame(n*Sa+1:row/alpha);zeros(L - (row/alpha-(n*Sa+1)+1),1)];
    end
    
    m = signal_in(n*Ss+1:L+(n*Ss));
    %window m
    windowed_m = w.*m;
    %find magnitude of m
    M = abs(fft(windowed_m));
    
    for i = 1:iterations  
        windowed_tmp = w.*temp_frame;
        RES = fft(windowed_tmp);
        phase = angle(RES) ;
        %recombine signal
        x = M.*exp(j*phase); 
        est_frame = ifft(x);
        temp_frame = real(est_frame);
    end
    est_frame = temp_frame;
    
    %update and combine partial frame for next iteration
    if (L+(n*Sa)<=row/alpha)
        par_frame(n*Sa+1:L+(n*Sa)) = par_frame(n*Sa+1:L+(n*Sa)) + est_frame;
    %make last frame in signal the right frame size 
    else
        par_frame(n*Sa+1:row/alpha) = par_frame(n*Sa+1:row/alpha) + est_frame(1:(row/alpha)-(n*Sa+1)+1);
    end
end
out = par_frame;
save output.mat out;
