function [x] = audioEcho(inVec,Fs,delay,echoGain)
%play audio clip
%delay the audio clip, lower the sound
%play the clip a second time with above changes

samples = round(Fs*delay);
ds = floor(samples);
signal = zeros(length(inVec)+ds,1);
signal(1:length(inVec))=inVec;
echo_signal = zeros(length(inVec)+ds,1);
echo_signal(ds+(1:length(inVec*echoGain))) = inVec*echoGain;
x = signal + echo_signal;
p = max(abs(x));
    if p>1
        x= x./ p;
    else
        x=x;
    end

end