function [BounceTime_T] = BounceTime(Trial)

%Counter
j = 1;
%This finds the minimum of the bounce arks
    for i = 3:length(Trial)-2
        if Trial(i,3) <= Trial(i+1,3) && Trial(i,3) <= Trial(i-1,3) && Trial(i,3) <= Trial(i+2,3) && Trial(i,3) <= Trial(i-2,3)

            BounceTime(j) = Trial(i,1);
            
            j=j+1;
        end
    end
% This makes sure that reference frames where the ball didnt quite hit the ball was found, a point
% was chosen in between where the ball is assumed to hit the ground.
    for z = 1:length(BounceTime)-1
       if round(BounceTime(z),1) == round(BounceTime(z+1),1)
          BounceTime(z) = (BounceTime(z) + BounceTime(z+1))/2;
          BounceTime(z+1) = [];
       elseif z == length(BounceTime)-1
           
           break
           
       end
    end
    %This gives us the actual time for each bounce instead of squential
    %time marks.
   for k = 1:length(BounceTime)-1
      BounceTime_T(k) = BounceTime(k+1)-BounceTime(k); 
       
   end
end