function fftBar(soundVec,Fs)
samples = 1000;
x = soundVec(1:samples); %first 1000 samples of the audio

% Make FFT
fftlength = 2^nextpow2(length(x));
ffft_one_sided_length = fftlength/2;
ffft = fft(x,fftlength);
ffft_normalised = ffft./(max(ffft));
ffft_shift = fftshift(ffft);
ffft_one_sided_normalised = ffft_normalised(1:ffft_one_sided_length);
fft_freq_x_axis_one_sided = Fs*(0:ffft_one_sided_length-1)/fftlength;

% plotting fft
figure
subplot(2,1,1)
plot(abs(ffft))
title('fft')
xlim([0 samples])

subplot(2,1,2)
plot(fft_freq_x_axis_one_sided./1000,abs(ffft_one_sided_normalised))
title('normalised and single-sided fft')
xlabel('frequency [KHz]')
ylabel('Amplitude')
xlim([0 5])

end