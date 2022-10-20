%Load data
S0 = load('key.mat');
S1 = load('studentID.mat');
S2 = load('studentAnswers.mat');
S3 = load('essayAnswers.mat');

%Loop for finding the scores from the multiple choice part of exam
for j=2:102

        for i=1:17
        LIA(i,1) = isequal(S0.key(1,i),S2.studentAnswers(j,i));
        end
 score(j,1)=sum(LIA);
end
MultScore = score(2:end,1);

weightedMultScore = MultScore(1:end,1).*(60/17);


%Finding the Essay portion of the exam
EssayScore = S3.essayAnswers(1:end,1);
weightedEssayScore = EssayScore .* (8/3);

Grade = cell(101,1);
Total = weightedMultScore + weightedEssayScore;

for i=1:101

if Total(i,1) > 93
        Grade{i,1} = 'A';
        elseif Total(i,1) > 90 
           Grade{i,1} = 'A-';
        elseif Total(i,1) >= 87 
           Grade{i,1} = 'B+';
        elseif Total(i,1) >= 83 
           Grade{i,1} = 'B';
        elseif Total(i,1) >= 80 
           Grade{i,1} = 'B-';
        elseif Total(i,1) >= 77 
           Grade{i,1} = 'C+';
        elseif Total(i,1) >= 73 
           Grade{i,1} = 'C';
        elseif Total(i,1) >= 70 
           Grade{i,1} = 'C-';
        elseif Total(i,1) >= 67 
           Grade{i,1} = 'D+';
        elseif Total(i,1) >= 63 
           Grade{i,1} = 'D';
        elseif Total(i,1) >= 60 
           Grade{i,1} = 'D-';
        elseif Total(i,1) < 60
       Grade{i,1} = 'F';
    
end
   
end

for i = 1:101
     
    fprintf(' Student %s got a %0.1f earning them the letter grade %s.\n',studentID{i,1},Total(i,1),Grade{i,1})  
   
end



    





 







