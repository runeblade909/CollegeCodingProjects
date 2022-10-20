function [x,v] = eigenmode4(L,ev,ne,nsub,scale)
%This function creates the eigenmode outputs from eigen values and
%frequencies
EigenValue4 = [0;0;0;0;ev]; %% For full 8x8
% ne=2; 
% nsub = 10; 
% scale = 1;
nv=ne*nsub+1; 
Le=L/ne; 
dx=Le/nsub; 
k=1;

x = zeros(nv,1);
v = zeros(nv,1); %% declare and set to zero plot arrays

for e = 1:ne  %% loop over elements
    
    xi = Le*(e-1); 
    vi = EigenValue4(2*e-1); 
    qi = EigenValue4(2*e); 
    vj = EigenValue4(2*e+1); 
    qj = EigenValue4(2*e+2);
    
    for n = 1:nsub %% loop over subdivisions
        xk = xi+dx*n; 
        z = (2*n-nsub)/nsub; %% isoP coordinate
        vk = (0.125*(4*(vi+vj)+2*(vi-vj)*(z^2-3)*z + Le*(z^2-1)*(qj-qi+(qi+qj)*z))); %% Hermitian interpolant
        k  = k+1;
        
        x(k)=xk;
        v(k)=-vk; %% build plot functions
    end  %% end n loop
    
end %% end e loop

v = v./mean(v);
v = scale*v;

end