close all;
clear all;
clc;

%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
%==========================================================================
%                  Magnetic tracking using a lookup table
%--------------------------------------------------------------------------
% Implemented for 4 coils
% You write the raw magnetic data for the five phases of the tracking in
% the "PARAMETERS" section and you can run the code to get the coordinates
% in cm of the sensor in the map that you loaded under the filename
% "lookup_table_filename"
% 
% Version: 1.0
% 
% Description:
% 1.0: start of the code
%==========================================================================
%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

%==========================================================================
%%                              PARAMETERS
%==========================================================================

lookup_table_filename="LUT_4_coils_spiral_v3_norm_already_calculated.txt";

%Enter your magnetic field data in mT
%raw values WITH the offset are expected as the offset is removed later on

% %used in the example: 241143,24.0,32.5,8.0,[0 0 0],[ 0.204731  0.036559 -0.004658],[ 0.014239  0.022884 -0.022161],[-0.434624  0.181093  0.507806],[-0.011188  0.041956 -0.028024]
% Bx_no_coil=0;
% By_no_coil=0;
% Bz_no_coil=0;
% 
% BxA_raw=0.204731;
% ByA_raw=0.036559;
% BzA_raw=-0.004658;
% 
% BxB_raw=0.014239;
% ByB_raw=0.022884;
% BzB_raw=-0.022161;
% 
% BxC_raw=-0.434624;
% ByC_raw=0.181093;
% BzC_raw=0.507806;
% 
% BxD_raw=-0.011188;
% ByD_raw=0.041956;
% BzD_raw=-0.028024;

%used in the example: 128005,12.5,36.5,16.5,[0 0 0],[0.030158 0.078411 0.145001],[ 0.00261   0.027671 -0.002096],[-0.054531  0.020254  0.017145],[-0.009429  0.014278 -0.004401]
% Bx_no_coil=0.051;
% By_no_coil=0;
% Bz_no_coil=0;
% 
% BxA_raw=-2.5;
% ByA_raw=-0.0;
% BzA_raw=-0.0;
% 
% BxB_raw=0.147;
% ByB_raw=-0.0;
% BzB_raw=-0.0;
% 
% BxC_raw=0.127;
% ByC_raw=0.000;
% BzC_raw=-0.0;
% 
% BxD_raw=0.077;
% ByD_raw=0.0;
% BzD_raw=-0.0;

Bx_no_coil=0.0378;
By_no_coil=-0.0541;
Bz_no_coil=-63.516;

BxA_raw=0.102;
ByA_raw=-0.0178;
BzA_raw=-137.679;

BxB_raw=-0.172;
ByB_raw=0.496;
BzB_raw=-137.586;

BxC_raw=0.109;
ByC_raw=-0.099;
BzC_raw=-139.216;

BxD_raw=-0.251;
ByD_raw=-0.688;
BzD_raw=-34.916;

%==========================================================================
%%                                MAIN
%==========================================================================
%% STEP 1
%Read the map
%The lookup table was changed so that the norm is already calculated and
%thus the code runs faster

[mapA,mapB,mapC,mapD,findXYZ,step_X,step_Y,step_Z,xmin,ymin,zmin] = load_mapping_LUT(lookup_table_filename);

%% STEP 2
%Prepare the maps for the tracking: we need to convert the 3D arrays to a
%1D array
%as reshape works only on 2D arrays, the 3D arrays are first convert to
%2D arrays

mapA_2D=zeros(length(mapA(:,1,1)),length(mapA(1,:,1))*length(mapA(1,1,:)));
mapB_2D=zeros(length(mapB(:,1,1)),length(mapB(1,:,1))*length(mapB(1,1,:)));
mapC_2D=zeros(length(mapC(:,1,1)),length(mapC(1,:,1))*length(mapC(1,1,:)));
mapD_2D=zeros(length(mapD(:,1,1)),length(mapD(1,:,1))*length(mapD(1,1,:)));

for i = 1:length(mapA_2D(:,1,1))
    for j = 1:length(mapA(1,:,1))
        mapA_2D(i,j*length(mapA(1,1,:))-length(mapA(1,1,:))+1:j*length(mapA(1,1,:)))=mapA(i,j,:);
        mapB_2D(i,j*length(mapB(1,1,:))-length(mapB(1,1,:))+1:j*length(mapB(1,1,:)))=mapB(i,j,:);
        mapC_2D(i,j*length(mapC(1,1,:))-length(mapC(1,1,:))+1:j*length(mapC(1,1,:)))=mapC(i,j,:);
        mapD_2D(i,j*length(mapD(1,1,:))-length(mapD(1,1,:))+1:j*length(mapD(1,1,:)))=mapD(i,j,:);
    end
end

%use of transposes as I want to flatten the array in a row-major order
%(C-style) (default in MATLAB is a column-major order (Fortran-style))
mapA_one_column=reshape(mapA_2D', [], 1)';
mapB_one_column=reshape(mapB_2D', [], 1)';
mapC_one_column=reshape(mapC_2D', [], 1)';
mapD_one_column=reshape(mapD_2D', [], 1)';
ones=ones(1,length(mapA_one_column));

%% STEP 3
%Prepare the magnetic field values

BxA=BxA_raw-Bx_no_coil;
ByA=ByA_raw-By_no_coil;
BzA=BzA_raw-Bz_no_coil;

BxB=BxB_raw-Bx_no_coil;
ByB=ByB_raw-By_no_coil;
BzB=BzB_raw-Bz_no_coil;

BxC=BxC_raw-Bx_no_coil;
ByC=ByC_raw-By_no_coil;
BzC=BzC_raw-Bz_no_coil;

BxD=BxD_raw-Bx_no_coil;
ByD=ByD_raw-By_no_coil;
BzD=BzD_raw-Bz_no_coil;

Btot=[[Bx_no_coil By_no_coil Bz_no_coil];
      [BxA        ByA        BzA];
      [BxB        ByB        BzB];
      [BxC        ByC        BzC];
      [BxD        ByD        BzD]];

%% STEP 4
%Launch the tracking and get the position x,y,z
[x, y, z] = calculs_of_points_LUT(Btot, 1, mapA_one_column, mapB_one_column, mapC_one_column, mapD_one_column, step_X, step_Y, step_Z, xmin, ymin, zmin, ones, findXYZ);

%results in cm
x,y,z