function [outVec, newFs] = compressMono(inVec, fs)
%Function for compressing audio vector
%by removing every other element

%Inputs:
%ratio = compression ratio
%inVec = input MONO audio vec
%fs = sampling frequency (Hz)
ratio = input('What ratio would you like for the compression?\n');


%create and indexing vector
indexVec = 1:ratio:length(inVec);

%What will be the length of the input vector
outVec = inVec(indexVec);

%Update the sampling freq.
newFs = fs/ratio;



end