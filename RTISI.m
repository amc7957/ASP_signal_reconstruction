function ser = RTISI(signal_in, fs, iter)
% signal_estimation(signal_in,fs,iter) implements the RTISI algorithm presented 
% in 'Real-Time Signal Estimation From Modified Short-Time Fourier
% Transform Magnitude Spectra' by Xinglei Zhu, Gerald Beauregard, and Lonce
% L. Wyse

[row,col] = size(signal_in);
%L window length
L = 1024;
%S step size 
S = L/4; 
a = 0.54;
b = -0.46;
n = 1:L;
w = (2*sqrt(S))/sqrt((4*a^2+2*b^2)*L)*(a+b*cos(2*pi*n/L)).';
iterations = iter;

par_frame = zeros(row,1);

for n = 0:((row-L)/S)
    
    %create partial frame
    if L+(n*S)<=row
        temp_frame = par_frame(n*S+1:L+(n*S));
    %make last frame in signal the right frame size    
    else
        temp_frame = [par_frame(n*S+1:row);zeros(L - (row-(n*S+1)+1),1)];
    end
    
    m = signal_in(n*S+1:L+(n*S));
    %window m
    windowed_m = w.*m;
    %take m into freq domain
    M = abs(fft(windowed_m));
    
    for i = 1:iterations
        windowed_tmp = w.*temp_frame;
        TEMP = fft(windowed_tmp);
        phase = angle(TEMP) ;
        %recombine signal
        x = M.*exp(j*phase); 
        est_frame = ifft(x);
        temp_frame = real(est_frame);
    end
    est_frame = temp_frame;
    
    %update and combine partial frame for next iteration
    if (L+(n*S)<=row)
        par_frame(n*S+1:L+(n*S)) = par_frame(n*S+1:L+(n*S)) + est_frame;
    %make last frame in signal the right frame size     
    else
        par_frame(n*S+1:row) = par_frame(n*S+1:row) + est_frame(1:row-(n*S+1)+1);
    end
end
out = par_frame;
save output.mat out;
ser = SER(signal_in(:,1),par_frame);
