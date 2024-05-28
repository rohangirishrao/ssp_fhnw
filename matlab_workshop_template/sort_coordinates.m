function [coordX_sorted, coordY_sorted, coordZ_sorted] = sort_coordinates(X, Y, Z)
    %For each axis, it gives only one examplary of each coordinate and
    %gives them in the increasing order

    coordX_sorted = unique(X);
    coordX_sorted = sort(coordX_sorted);
    
    coordY_sorted = unique(Y);
    coordY_sorted = sort(coordY_sorted);
    
    coordZ_sorted = unique(Z);
    coordZ_sorted = sort(coordZ_sorted);
end