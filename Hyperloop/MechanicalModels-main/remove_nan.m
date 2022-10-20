function new_data = remove_nan(data)

    [r, c] = size(data);
    new_data = zeros(r,c);
    
    for i = 1:r
        for j = 1:c
            
            
            
            if isnan(data(i,j))
                new_data(i,j) = (data(i-1,j) + data(i+1,j))/2;
            else
                new_data(i,j) = data(i,j);
            end
        end
    end
end

