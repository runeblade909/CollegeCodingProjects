


function[Q,t] = easterSunday(y)
%y/19, remainder =a. 
% y/100, quotient=b and remainder=c
% b/4, quotient=d and remainder=e.
% (8*b+13)/25, quotient=g
%(19*a+b-d-g+15)/30, remainder=h
% c/4, quotient=j and remainder=k
% (a+11*h)/319, quotient=m
% (2*e+2*j-k-h+m+32)/7, remainder=r
% (h-m+r+90)/25, quotient=n
% (h-m+r+n+19)/32, remainder=p
[q,r] = quorem(sym(y),sym(19));
out=[q,r];
end