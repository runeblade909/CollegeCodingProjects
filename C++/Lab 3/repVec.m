function [returnVec] = repVec(x,y)
clc
%Hey, I did this in class when it was assigned, but I guess it didnt
%successfully submit when I submitted it. I understand its extremely late
%but I just looked and saw it was missing. If I could get any credit I
%would be appreciative.
returnVec=[];
 a=1;
 b=1;
  if isvector(x)==0
      error('I need a vector!')
  elseif y<0
      error('Negative number dont work chief')
      
  elseif iscolumn(x)==1
      ctranspose(x);
  end
      
for i= 1:length(x);
    for j =1:y;
        returnVec(b)=x(a);
        b=b+1;
    end
    a=a+1;
end



end