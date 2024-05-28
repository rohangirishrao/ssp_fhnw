function [x, y, z] = calculs_of_points_LUT(Btot, number_of_sensors, mapA_one_column, mapB_one_column, mapC_one_column, mapD_one_column, step_X, step_Y, step_Z, xmin, ymin, zmin, ones, findXYZ)
    %the tracking is done on the norm of the magnetic field, hence the
    %calculus of NormQuadruplet
    NormQuadruplet = zeros(number_of_sensors, 4);
    for k = 1:number_of_sensors
        NormQuadruplet(k, 1) = norm(Btot(2, (3*k-2):(3*k)));
        NormQuadruplet(k, 2) = norm(Btot(3, (3*k-2):(3*k)));
        NormQuadruplet(k, 3) = norm(Btot(4, (3*k-2):(3*k)));
        NormQuadruplet(k, 4) = norm(Btot(5, (3*k-2):(3*k)));
    end
    
    x = zeros(1,number_of_sensors);
    y = zeros(1,number_of_sensors);
    z = zeros(1,number_of_sensors);
    
    for k = 1:number_of_sensors
        position = find_coordinates_LUT(NormQuadruplet(k,:), mapA_one_column, mapB_one_column, mapC_one_column, mapD_one_column, step_X, step_Y, step_Z, xmin, ymin, zmin, ones, findXYZ);
        x(1,k)=position(1);
        y(1,k)=position(2);
        z(1,k)=position(3);
    end
end