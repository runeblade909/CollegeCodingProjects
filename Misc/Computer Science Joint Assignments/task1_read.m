clc
clear

%specify the file
fileName ='quiz_results.xlsx'; 


[nums,txt,raw] = xlsread(fileName);
%clean the data

key = txt(2,7:end);
studentID = txt(3:end,2);
studentAnswers = txt(2:end,7:end);
essayAnswers = nums(2:end,3);
%save data

save('key.mat','key');
save('studentID.mat','studentID');
save('studentAnswers.mat','studentAnswers');
save('essayAnswers.mat','essayAnswers');