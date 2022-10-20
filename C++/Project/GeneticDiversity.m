clc;clear;close all;
fileName = 'Fitness.csv';
fileName1 = 'max_ele.csv';
[~,~,RAW] = xlsread(fileName);
[~,~,MAX] = xlsread(fileName1);
[Generations,elements] = size(MAX);
[~,elements1] = size(RAW);
MAXARR = cell2mat(MAX);
RAWARR = cell2mat(RAW);
xValues = 1:Generations;
yValues = 0:1;



Avg = zeros(Generations,1);

% for j = 1:Generations
%     
%        Avg(j,1) = sum(RAW,j,elements1);
% 
%  
% end

Avg = ((sum(RAWARR,2))/elements1);


plot(xValues,Avg,xValues,MAXARR);
title('Genetic Diversity');
xlabel('Generations');
ylabel('Fitness');
xline(Generations);
legend({'Average','Max','Last Generation'}, 'Location','southeast');
