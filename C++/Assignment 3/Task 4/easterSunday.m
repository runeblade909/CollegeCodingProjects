%By Zak Reichenbach
%9/19/2019
function[Q,R] = easterSunday(y)
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
y = input('Year:');
[~,a] = qrFunc(y,19);
[b,c] = qrFunc(y,100);
[d,e] = qrFunc(b,4);
[g,~] = qrFunc(((8*b)+13),25);
[~,h] = qrFunc((19*a+b-d-g+15), 30);
[j,k] = qrFunc(c,4);
[m,~] = qrFunc((a+11*h),319);
[~,r] = qrFunc((2*e+2*j-k-h+m+32),7);
[n,~] = qrFunc((h-m+r+90),25);
[~,p] = qrFunc((h-m+r+n+19),32);
fprintf("In %d, Easter Sunday is on %d/%d.\n", y,n,p);

end
