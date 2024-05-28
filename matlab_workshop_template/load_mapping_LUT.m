function [mapA,mapB,mapC,mapD,findXYZ,step_X,step_Y,step_Z,xmin,ymin,zmin] = load_mapping_LUT(filename)
    %get number of lines in the file
    txt_file=fopen("LUT_4_coils_spiral_v3_norm_already_calculated.txt",'r');
    number_of_lines = 0;
    while ~feof(txt_file)
        line = fgetl(txt_file);
        if ~isempty(line)
            number_of_lines = number_of_lines + 1;
        end
    end
    fclose(txt_file);

    %load the file according to some parameters (types, delimiter etc.)
    opts = delimitedTextImportOptions("NumVariables", 9);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double"];

    data=readtable("LUT_4_coils_spiral_v3_norm_already_calculated.txt",opts);
    data=table2array(data);

    clear opts

    new_coordX=zeros(number_of_lines-1,1);
    new_coordY=zeros(number_of_lines-1,1);
    new_coordZ=zeros(number_of_lines-1,1);
    Ba_all_mappings=zeros(number_of_lines-1,1);
    Bb_all_mappings=zeros(number_of_lines-1,1);
    Bc_all_mappings=zeros(number_of_lines-1,1);
    Bd_all_mappings=zeros(number_of_lines-1,1);

    new_coordX(:,1)=data(:,2);
    new_coordY(:,1)=data(:,3);
    new_coordZ(:,1)=data(:,4);
    Ba_all_mappings(:,1)=data(:,6);
    Bb_all_mappings(:,1)=data(:,7);
    Bc_all_mappings(:,1)=data(:,8);
    Bd_all_mappings(:,1)=data(:,9);
    
    %The aim of the function "sort_coordinates" is to get all the
    %coordinates inside the LUT but only in one examplary and in the
    %increasing order
    [coordX_sorted,coordY_sorted,coordZ_sorted]=sort_coordinates(new_coordX,new_coordY,new_coordZ);

    mapA = zeros(length(coordX_sorted), length(coordY_sorted), length(coordZ_sorted));
    mapB = zeros(length(coordX_sorted), length(coordY_sorted), length(coordZ_sorted));
    mapC = zeros(length(coordX_sorted), length(coordY_sorted), length(coordZ_sorted));
    mapD = zeros(length(coordX_sorted), length(coordY_sorted), length(coordZ_sorted));

    findXYZ=zeros(number_of_lines-1,3);
    for x = 1:length(coordX_sorted)
        for y = 1:length(coordY_sorted)
            for z = 1:length(coordZ_sorted)
                ind=z+length(coordZ_sorted)*(y-1)+length(coordZ_sorted)*length(coordY_sorted)*(x-1);
                
                %will be used later to come back to the indices x, y
                %and z that gave the index "ind"
                findXYZ(ind,:)=[x,y,z];

                %fill the LUT in matrices
                mapA(x,y,z) = Ba_all_mappings(ind,1);
                mapB(x,y,z) = Bb_all_mappings(ind,1);
                mapC(x,y,z) = Bc_all_mappings(ind,1);
                mapD(x,y,z) = Bd_all_mappings(ind,1);
            end
        end
    end

    %calculate some parameters from the LUT data that will allow later to
    %convert the indices x, y and z found thanks to findXYZ to coordinates
    %X, Y and Z in cm
    xmin = min(new_coordX);
    ymin = min(new_coordY);
    zmin = min(new_coordZ);

    xmax = max(new_coordX);
    ymax = max(new_coordY);
    zmax = max(new_coordZ);

    step_X=(xmax-xmin)/(length(coordX_sorted)-1);
    step_Y=(ymax-ymin)/(length(coordY_sorted)-1);
    step_Z=(zmax-zmin)/(length(coordZ_sorted)-1);
end