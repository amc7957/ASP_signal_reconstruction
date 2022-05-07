function TSM(signal_in, alpha)
% TSM(signal_in, alpha) stretches a signal 'signal in' by a factor of
% alpha. This method is based on the OLA-TSM method presented in 'A Review
% of Time-Scale Modification of Music Signals' by Jonathan Dreidger


signal_in = signal_in.';

N = 5000;
Hs= N/2;
Ha = Hs/alpha;

[row,col] = size(signal_in);

shift_matrix = zeros(N,col);

for m = 0:((col-N)/Ha)-1
    %separate into analysis frames
    xm = signal_in(:,(m*Ha)+1:N+(m*Ha));
    %window
    w = hann(N).';
    ym = xm.*w;
    %%shift frame based on Hs
    shift_matrix(m+1,m*Hs+1:(m*Hs)+N) = ym;
end

%reconstruct signal
y_final = sum(shift_matrix,1);
save output.mat y_final;



