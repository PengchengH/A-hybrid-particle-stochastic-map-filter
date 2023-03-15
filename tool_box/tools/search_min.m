function [min_address,sed_min_address] = search_min(invec)

    if ~size(invec,2)>1
        error('Input should be a vector');
    end

    [value,address] = sort(invec);
    
    for k = 2:size(value,2)
        if value(k) > value(1)
            break;
        end
    end
    
    min_address = address(1:k-1);

    if size(min_address,2) ==1
        % if only one minimun, search second minmum
        for k = 3:size(value,2)
            if value(k) > value(2)
                break;
            end
        end 
        sed_min_address = address(2:k-1);

    elseif size(min_address,2) >1
        % if there are more than one minimum, no second minimum
        sed_min_address = zeros(0,0);         
    end

    

end

