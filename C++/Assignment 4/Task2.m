%Task2.m
%By: Zak Reichenbach
%9/26/2019

clc
clear

V = [11,5,3,2,-18,4,-5,5,-66]
%use find() function to get rid of negative values
X=find(V>0);

V1=V(X)
%use logicalVec to get rid of negative values
logicalVec= V<0;

V2=V(logicalVec==0)
