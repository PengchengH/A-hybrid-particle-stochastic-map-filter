function [num_vry_index] = parallel_arrange(num_base)

% variable number set
num_vry  = ones(1,size(num_base,2));
num_vry(size(num_base,2)) = 0;

i_sum    = prod(num_base);

num_vry_index = zeros(i_sum,size(num_base,2));

for i_index = 1:i_sum
    num_vry(size(num_base,2)) = num_vry(size(num_base,2))+1;
    for i = 1:size(num_base,2)
        j = size(num_base,2)+1-i;
        if num_vry(j) > num_base(j)
           num_vry(j)   = 1;
           num_vry(j-1) = num_vry(j-1)+1;
        end
    end
    num_vry_index(i_index,:) = num_vry;
end

end

