function [Bounce_T] = BounceHeight(Trial)

%Counter
j = 1;
for i = 3:length(Trial)-2
    %This checks for max values on graphs.
    if Trial(i,3) > Trial(i-2,3) && Trial(i,3) > Trial(i+2,3)
       Bounce(j) = Trial(i,3);
        j=j+1;
        k=i;
    elseif Trial(i,3) == Trial(i+1,3) 
     Bounce(j) = Trial(i,3);
        j=j+1;
        k=i;
    elseif i == length(Trial)-2
    break
    end
end

%This next loop checks for false max values on the loop;
%Counter
q=1;
for i = 1:length(Bounce)-1
   if round(Bounce(i),1) == round(Bounce(i+1),1) && Bounce(i) > Bounce(i+1)
       Bounce_T(q) = Bounce(i);
       q=q+1;
   elseif round(Bounce(i),1) == round(Bounce(i+1),1) && Bounce(i) < Bounce(i+1)
      Bounce_T(q) =Bounce(i+1); 
       q=q+1;
   end
end
%This further checks for duplicate max numbers
for i = 2:length(Bounce_T)
       if Bounce_T(i) == Bounce_T(i-1)
           Bounce_T(i) = [];
       elseif i == length(Bounce_T)
           break
       end
end
end
