rectangle = load("11_in_rectangle.mat");
sqaure = load("11inch_SquareRod.mat");


[m,n] = size(data);

onefive =  0;
oneeight  = 0 ;
onetwo = 0;


for i = 1:m

    if data(i,3) == 1.5
        onefive = onefive +  1;


    elseif data(i,3) == 1.875
        oneeight = oneeight + 1;

    elseif data(i,3) ==  1.25

        onetwo = onetwo + 1;

    end



end




idx = 1;

while data(idx,3) ~= 1.5

    idx = idx+1;
end

idxEnd = idx + 23;

data(idx:idxEnd,3) = 1.375;




onefive =  0;
oneeight  = 0 ;
onetwo = 0;


for i = 1:m

    if data(i,3) == 1.5
        onefive = onefive +  1;


    elseif data(i,3) == 1.875
        oneeight = oneeight + 1;

    elseif data(i,3) ==  1.25

        onetwo = onetwo + 1;

    end



end

