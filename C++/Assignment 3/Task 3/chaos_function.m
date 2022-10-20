%chaos_function.m
function[out]= chaos_function(start_number,step_number,end_number,vector)

vector=(start_number:step_number:end_number)
a=length(vector);
out=vector(a-1)-vector(2);
end
