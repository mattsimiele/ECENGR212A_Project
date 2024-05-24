function [num_out, den_out, transMatrix] = filter_z_transform(filter_type, num_in, den_in, c1, c2)
    % Computes the numerator, denominator, and transformation matrix for a given 
    % digital filter transformation 
    
    % Check if num_in and den_in have the same order
    if length(num_in) ~= length(den_in)
        error('Transfer function does not have the same order in den and numerator!');
    else 
        order = length(num_in)-1;
    end

    if strcmp(filter_type, 'lp')
        [num_trans, den_trans] = lp2lp_z_transform(c1);
    elseif strcmp(filter_type, 'hp')
        [num_trans, den_trans] = lp2hp_z_transform(c1);
    elseif strcmp(filter_type, 'bp')
        [num_trans, den_trans] = lp2bp_z_transform(c1, c2);
    elseif strcmp(filter_type, 'bs')
        [num_trans, den_trans] = lp2bs_z_transform(c1, c2);
    else
        error('Desired filter type is not supported. Please use: lp, hp, bp, or bs!');
    end

    transMatrix = compute_transfer_matrix(num_trans, den_trans, order);
    
    [num_out, den_out] = compute_output_transfer_func(transMatrix, num_in, den_in);
end

function [num_trans, den_trans] = lp2lp_z_transform(c1)
    num_trans = [-1*c1, 1];
    den_trans = [1, -1*c1];
end

function [num_trans, den_trans] = lp2hp_z_transform(c1)
    num_trans = [-1*c1, -1];
    den_trans = [1, c1];
end

function [num_trans, den_trans] = lp2bp_z_transform(c1, c2)
    num_trans = -1.*[(c2-1), -2*c1*c2, (c2+1)];
    den_trans = [(c2+1), -2*c1*c2, (c2-1)];
end

function [num_trans, den_trans] = lp2bs_z_transform(c1, c2)
    num_trans = [(1-c2), -2*c1, (1+c2)];
    den_trans = [(1+c2), -2*c1, (1-c2)];
end

function transMatrix = compute_transfer_matrix(num, den, order)
    transMatrix = zeros(order+1, (order * (length(num) - 1)) + 1);
    for i = 1:order+1
        num_conv = 1;
        den_conv = 1;
        for j = 1:i-1
            num_conv = conv(num_conv, num);
        end
        for j = 1:order - i + 1
            den_conv = conv(den_conv, den);
        end
        transMatrix(i,:) = conv(num_conv, den_conv);
    end
end

function [num_out, den_out] = compute_output_transfer_func(transMatrix, num, den)
    num_out = zeros(size(transMatrix(1,:)));
    den_out = zeros(size(transMatrix(1,:)));
    for i = 1:length(num)
       den_out = den_out + den(i) * transMatrix(i,:);
       num_out = num_out + num(i) * transMatrix(i,:);  
    end
    num_out = num_out / den_out(1);
    den_out = den_out / den_out(1);
end
