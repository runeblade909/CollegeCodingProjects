clc
clear


X=randi([0 5],4,5)


y=randi(20)

t = logical(X)==0

newVec= X(t)

z=length(newVec)

if y > z
    q = X > 3

yeeVec=X(q)

e=sum(yeeVec)
elseif y <= z
    q =max(X)

end
    
