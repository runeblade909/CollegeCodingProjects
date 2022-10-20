%audioCompress.m
clear
clc

%audioread() function for getting audio into
%vector (or matrix) format

[y,fs] = audioread('talkieMono.wav');
%sound(y,fs)

compressionRatio = input(' Enter the desired compression ratio:');

[compressedVec, fs] = compressMono(compressionRatio,y,fs);

sound(compressedVec,fs)

audiowrite('compressed.wav',compressedVec,fs);