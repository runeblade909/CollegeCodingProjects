function [C] = makemat(A,B)
clc

if isrow(A)
    transpose(A);
elseif isrow(B)
    transpose(B);
end

if length(A)==length(B)
    X = A ;
    Y = B ;
elseif length(B)>length(A)
  X = zeros(size(B)) ;
  X(1:length(A)) = A ;
  Y = B ;
elseif length(A)>length(B)
    Y = zeros(size(A)) ;
    Y(1:length(B)) = B ;
    X = A ;
end
C(:,1)=X ;
C(:,2)=Y ;




end