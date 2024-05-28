function [deduced_position_XYZ] = find_coordinates_LUT(Bmeasure,mapA_one_column, mapB_one_column, mapC_one_column, mapD_one_column, step_X, step_Y, step_Z, xmin, ymin, zmin,ones,findXYZ)
    %Calculation of the cost function (without the norm)
    all_values = [ones - mapA_one_column./Bmeasure(1), ones - mapB_one_column./Bmeasure(2), ones - mapC_one_column./Bmeasure(3), ones - mapD_one_column./Bmeasure(4)];
    
    %reshape to have the values concerning a coordinate on the same
    %line
    tmp_to_calculate_norms=reshape(all_values, [], 4);
    %calculate the norm line per line
    norms=vecnorm(tmp_to_calculate_norms,2,2);
    
    %find the index of the minimum of the cost function
    [~, ind] = min(norms);
    
    %get the indices of the minimum in the LUT
    x = findXYZ(ind, 1);
    y = findXYZ(ind, 2);
    z = findXYZ(ind, 3);
    
    %convert these indices to real coordinates in cm
    deduced_position_XYZ = [(x-1)*step_X+xmin, (y-1)*step_Y+ymin, (z-1)*step_Z+zmin];
end