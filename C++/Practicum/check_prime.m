function[is_prime] = check_prime (a_num)
    clc

    if isscalar(a_num) ~= 1
        fprintf('Thats not a number bub.\n')
    end
is_prime = 1;
%start loop: iterator 'i' takes values from 2 till (a_num - 1). If the
%modulus(or remainder) of a_num with 'i' is equal to 0, then set is_prime=0

i=2:(a_num - 1);
x=rem(a_num,i);
y = x == 0;



    if sum(y) == 0
        fprintf (' This be a prime number.\n')
        is_prime = 1;
    else 
        fprintf (' This aint no prime number.\n')
        is_prime = 0;
    end

end