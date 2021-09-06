function generate_cable_data(filename, length, fmin, npc, npd)
% filename : name of data file to generate
% length : cable length in meters
% fmin : minimal frequency for computation
% npd : Number of points per decade for frequency scale
% npc : Number of decades for frequency scale

    [~,filename_noext] = fileparts(filename);
    mkdir('femm');
    mkdir('results');
    
    fmin_decade = fix(log(fmin)/log(10));
    freq = logspace(fmin_decade, fmin_decade+npc, npd * npc);
    
    freq = fmin;
    for i=2:(npd*npc+1)
        freq(i) = freq(i-1) * 10^(1/npd);
    end
    
    [Lf, Rf, C, Rdc, nphs] = femm_computation(length, freq);
    %load(strcat('.\results\',filename_noext)); % for debug purpose
    %nphs = 4;
    % per-unit computation
    Lf_km = Lf / (length/1000);
    Rf_km = Rf / (length/1000);
    C_km = C / (length/1000);
    Rdc_km = Rdc / (length/1000);
    
    save(strcat('.\results\',filename_noext), 'Lf', 'Rf', 'C', 'length', 'fmin', 'npd', 'npc', 'freq', 'Rdc');
    
    %% Write cable model data
    fileID = fopen(filename,'w');
%     fprintf(fileID,'C        1         2         3         4         5         6         7         8\n');
%     fprintf(fileID,'C 345678901234567890123456789012345678901234567890123456789012345678901234567890\n');
%     fprintf(fileID,'BEGIN NEW DATA CASE\nCABLE CONSTANTS\n');
%     fprintf(fileID,'CABLE-MODEL         FD-MODEL  FDQ             1000%10g\n',length);
%     fprintf(fileID,'EXTERNAL FILENAME\n');
%     fprintf(fileID,'EXT_CABLE_DATA.DAT\n');
%     fprintf(fileID,'%10d%10g%10g%10d%10d\n', nphs, 0.001, fmin, npd, npc);
%     fclose(fileID);
    
    %% Write external cable data file
    for i= 1:1:max(size(freq))
        f = freq(1,i);
        Y(:,:,i) = admittance(f,C_km); % Calculation for Y(admittance) matrix
        Z(:,:,i) = impedance(f,Lf_km(:,:,i),Rf_km(:,:,i)); % Calculation for Z (impedence) matrix
    end
    
    %fileID = fopen('EXT_CABLE_DATA.DAT','w');
    fprintf(fileID,'C        1         2         3         4         5         6         7         8\n');
    fprintf(fileID,'C 345678901234567890123456789012345678901234567890123456789012345678901234567890\n');
    fprintf(fileID,'C Number of phases\n%d\n', nphs);
    
    lines_number = fix(nphs * 2 / 5);
    columns_remaining = nphs - lines_number;
    
    i=1;
    fprintf(fileID,'%G',freq(i));
    idx = 1;
    size_mat = size(Y(:,:,1));
    for j = 1:size_mat(1)
        for k = 1:size_mat(2)
            if (j>=k)
            Y_str(idx) = real(Y(j,k,i));
            Y_str(idx+1) = imag(Y(j,k,i));
            Z_str(idx) = real(Z(j,k,i));
            Z_str(idx+1) = imag(Z(j,k,i));
            idx = idx + 2;
            end
        end
    end

    for j = 1:max(size(Z_str))
        if (mod(j,5) == 1)
            fprintf(fileID,'\n');
        else
            fprintf(fileID,' ');
        end
        fprintf(fileID, '%G', Z_str(j));
    end
    %fprintf(fileID,'\n');

    for j = 1:max(size(Y_str))
        if (mod(j,5) == 1)
            fprintf(fileID,'\n');
        else
            fprintf(fileID,' ');
        end
        fprintf(fileID, '%G', Y_str(j));

    end
    fprintf(fileID,'\n');
    
    fprintf(fileID, '0');
    idx = 1;
    size_mat = size(Y(:,:,1));
    for j = 1:size_mat(1)
        for k = 1:size_mat(2)
            if (j>=k)
            Y_str(idx) = 0;
            Y_str(idx+1) = 0;
            Z_str(idx) = Rdc_km(j,k);
            Z_str(idx+1) = 0;
            idx = idx + 2;
            end
        end
    end

    for j = 1:max(size(Z_str))
        if (mod(j,5) == 1)
            fprintf(fileID,'\n');
        else
            fprintf(fileID,' ');
        end
        fprintf(fileID, '%G', Z_str(j));

    end
    %fprintf(fileID,'\n');

    for j = 1:max(size(Y_str))
        if (mod(j,5) == 1)
            fprintf(fileID,'\n');
        else
            fprintf(fileID,' ');
        end
        fprintf(fileID, '%G', Y_str(j));
    end
    fprintf(fileID,'\n');
    
    for i = 1:max(size(freq))
        fprintf(fileID,'%G',freq(i));
        idx = 1;
        size_mat = size(Y(:,:,1));
        for j = 1:size_mat(1)
            for k = 1:size_mat(2)
                if (j>=k)
                Y_str(idx) = real(Y(j,k,i));
                Y_str(idx+1) = imag(Y(j,k,i));
                Z_str(idx) = real(Z(j,k,i));
                Z_str(idx+1) = imag(Z(j,k,i));
                idx = idx + 2;
                end
            end
        end
        
        for j = 1:max(size(Z_str))
            if (mod(j,5) == 1)
                fprintf(fileID,'\n');
            else
                fprintf(fileID,' ');
            end
            fprintf(fileID, '%g', Z_str(j));
        end
        %fprintf(fileID,'\n');
        
        for j = 1:max(size(Y_str))
            if (mod(j,5) == 1)
                fprintf(fileID,'\n');
            else
                fprintf(fileID,' ');
            end
            fprintf(fileID, '%g', Y_str(j));
        end
        fprintf(fileID,'\n');
        
    end
    fclose(fileID);
   
end

function [Lf, Rf, C, Rdc, nphs] = femm_computation(length, freq)
%%% matlab code for parameter calculation L, R and C with FEMM %%%
    nphs = 4;
    l=length; %per unit km calculation
    %len= 100; %length of cable
    % freq =[1.00000e-05, 1.00000e-02, 1.25893e-02, 1.58489e-02, 1.99526e-02, 2.51189e-02, 3.16228e-02, 3.98107e-02, 5.01187e-02, 6.30957e-02, 7.94328e-02, 1.00000e-01, 1.25893e-01, 1.58489e-01, 1.99526e-01, 2.51189e-01, 3.16228e-01, 3.98107e-01, 5.01187e-01, 6.30957e-01, 7.94328e-01, 1.00000e+00, 1.25893e+00, 1.58489e+00, 1.99526e+00, 2.51189e+00 , 3.16228e+00, 3.98107e+00, 5.01187e+00, 6.30957e+00, 7.94328e+00, 1.00000e+01, 1.25893e+01, 1.58489e+01, 1.99526e+01, 2.51189e+01, 3.16228e+01, 3.98107e+01, 5.01187e+01, 6.30957e+01, 7.94328e+01, 1.00000e+02, 1.25893e+02, 1.58489e+02, 1.99526e+02, 2.51189e+02, 3.16228e+02, 3.98107e+02, 5.01187e+02, 6.30957e+02, 7.94328e+02, 1.00000e+03, 1.25893e+03, 1.58489e+03, 1.99526e+03, 2.51189e+03, 3.16228e+03, 3.98107e+03, 5.01187e+03, 6.30957e+03, 7.94328e+03, 1.00000e+04, 1.25893e+04, 1.58489e+04, 1.99526e+04, 2.51189e+04, 3.16228e+04, 3.98107e+04, 5.01187e+04, 6.30957e+04, 7.94328e+04, 1.00000e+05, 1.25893e+05, 1.58489e+05, 1.99526e+05,  2.51189e+05, 3.16228e+05, 3.98107e+05, 5.01187e+05, 6.30957e+05, 7.94328e+05, 1.00000e+06];
    %  freq =[1.00000e-05, 1.00000e-02, 1.25893e-02, 1.58489e-02, 1.99526e-02];
    %   freq= [logspace(-5,0,20) logspace(0,2,20) logspace(2,3,100) logspace(3,4,200)];
    %freq=50;
 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Calculation of Rdc matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rho1= 1.728e-8;
    r1= 0.032; % meter
    s1= pi*(r1)^2; %metersqaure
    R11= (rho1*l)/s1;
    R33= R11;

    rho2= 2.83e-8;
    ri2= 0.0569; % meter
    ro2= 0.0582;
    si2= pi*(ri2)^2; %metersqaure
    so2= pi*(ro2)^2; %metersqaure
    s2 = so2-si2;
    R22= (rho2*l)/s2;
    R44= R22;

    Rdc = diag([R11,R22,R33,R44]);
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       SIMULATION for C matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    openfemm
    newdocument(1);

    %Définition du pb : 
    ei_probdef('meters','planar',1.e-8,l,30);

    ei_getmaterial('Air');
    ei_addmaterial('soil', 3, 3, 0);
    ei_addmaterial('insulation', 2.5, 2.5, 0); % ε permitivity= 2.5

    %cable geometry-1
    %core
    ei_addnode(0.032,0);
    ei_addnode(-0.032,0);
    ei_addarc(0.032,0,-0.032,0,180,1);
    ei_addarc(-0.032,0,0.032,0,180,1);
    ei_addblocklabel(0,0); 
    ei_selectlabel(0,0);
    ei_setblockprop('Air', 0, 0, 0);
    ei_clearselected;

    ei_selectarcsegment(0,0.01);
    ei_setarcsegmentprop(1,0, 0, 0, 'v1');
    ei_clearselected;
    ei_selectarcsegment(0,-0.01);
    ei_setarcsegmentprop(1,0, 0, 0, 'v1');
    ei_clearselected;

    %insulalation
    ei_addnode(0.0569,0);
    ei_addnode(-0.0569,0);
    ei_addarc(0.0569,0,-0.0569,0,180,1);
    ei_addarc(-0.0569,0,0.0569,0,180,1);
    ei_addblocklabel(-0.04445,0);
    ei_selectlabel(-0.04445,0);
    ei_setblockprop('insulation', 0, 0, 0);
    ei_clearselected;

    %sheath
    ei_addnode(0.0582,0);
    ei_addnode(-0.0582,0);
    ei_addarc(0.0582,0,-0.0582,0,180,1);
    ei_addarc(-0.0582,0,0.0582,0,180,1);
    ei_addblocklabel(0.05755,0); 
    ei_selectlabel(0.05755,0);
    ei_setblockprop('Air', 0, 0, 0);
    ei_clearselected;

    ei_selectarcsegment(0,0.05);
    ei_setarcsegmentprop(1,0, 0, 0, 'v2');
    ei_clearselected;
    ei_selectarcsegment(0,-0.05);
    ei_setarcsegmentprop(1,0, 0, 0, 'v2');
    ei_clearselected;

    % outer insulalation
    ei_addnode(0.0639,0);
    ei_addnode(-0.0639,0);
    ei_addarc(0.0639,0,-0.0639,0,180,1);
    ei_addarc(-0.0639,0,0.0639,0,180,1);
    ei_addblocklabel(-0.06105,0);
    ei_selectlabel(-0.06105,0);
    ei_setblockprop('insulation', 0, 0, 0);
    ei_clearselected;

    %cable geometry-2
    %core
    ei_addnode(0.468,0);
    ei_addnode(0.532,0);
    ei_addarc(0.468,0,0.532,0,180,1);
    ei_addarc(0.532,0,0.468,0,180,1);
    ei_addblocklabel(0.5,0); 
    ei_selectlabel(0.5,0);
    ei_setblockprop('Air', 0, 0, 0);
    ei_clearselected;

    ei_selectarcsegment(0.5,0.01);
    ei_setarcsegmentprop(1,0, 0, 0, 'v3');
    ei_clearselected;
    ei_selectarcsegment(0.5,-0.01);
    ei_setarcsegmentprop(1,0, 0, 0, 'v3');
    ei_clearselected;

    %insulalation
    ei_addnode(0.4431,0);
    ei_addnode(0.5569,0);
    ei_addarc(0.4431,0,0.5569,0,180,1);
    ei_addarc(0.5569,0,0.4431,0,180,1);
    ei_addblocklabel(0.45555,0);
    ei_selectlabel(0.45555,0);
    ei_setblockprop('insulation', 0, 0, 0);
    ei_clearselected;

    %sheath
    ei_addnode(0.5582,0);
    ei_addnode(0.4418,0);
    ei_addarc(0.5582,0,0.4418,0,180,1);
    ei_addarc(0.4418,0,0.5582,0,180,1);
    ei_addblocklabel(0.55755,0); 
    ei_selectlabel(0.55755,0);
    ei_setblockprop('Air', 0, 0, 0);
    ei_clearselected;

    ei_selectarcsegment(0.5,0.05);
    ei_setarcsegmentprop(1,0, 0, 0, 'v4');
    ei_clearselected;
    ei_selectarcsegment(0.5,-0.05);
    ei_setarcsegmentprop(1,0, 0, 0, 'v4');
    ei_clearselected;

    % outer insulalation
    ei_addnode(0.5639,0);
    ei_addnode(0.4361,0);
    ei_addarc(0.5639,0,0.4361,0,180,1);
    ei_addarc(0.4361,0,0.5639,0,180,1);
    ei_addblocklabel(0.43895,0);
    ei_selectlabel(0.43895,0);
    ei_setblockprop('insulation', 0, 0, 0);
    ei_clearselected;

    %définition des conductors
    ei_addconductorprop('v1', 1, 0, 1);
    ei_addconductorprop('v2', 0, 0, 1);
    ei_addconductorprop('v3', 0, 0, 1);
    ei_addconductorprop('v4', 0, 0, 1);

    %soil line
    ei_drawline(-50,1.5,50.5,1.5);
    ei_addblocklabel(0.25,1);
    ei_selectlabel(0.25,1);
    ei_setblockprop('soil', 0, 0, 0);
    ei_clearselected;

    %définition des conditions aux limites
    ei_addboundprop('boundary', 0, 0, 0, 0, 0);
    ei_drawline(-50,50,-50,-50);
    ei_drawline(50.5,50,50.5,-50);
    ei_drawline(-50,-50,50.5,-50);
    ei_drawline(-50,50,50.5,50)

    ei_selectsegment(-50,48);
    ei_setsegmentprop('boundary',0,0,0,0,0);
    ei_clearselected;
    ei_selectsegment(-48,-50);
    ei_setsegmentprop('boundary',0,0,0,0,0);
    ei_clearselected;
    ei_selectsegment(50.5,-48);
    ei_setsegmentprop('boundary',0,0,0,0,0);
    ei_clearselected;
    ei_selectsegment(-48,50);
    ei_setsegmentprop('boundary',0,0,0,0,0);
    ei_clearselected;

    ei_addblocklabel(0.25,8);
    ei_selectlabel(0.25,8);
    ei_setblockprop('Air', 0, 0, 0);
    ei_clearselected;

    ei_zoomnatural
    ei_saveas('.\femm\cableY.FEE')

    %Results of simulation

    ei_analyze(0)
    ei_loadsolution

    values1= eo_getconductorproperties('v1');
    C11= values1(1,2)/values1(1,1);
    values12 = eo_getconductorproperties('v2');
    C12= values12(1,2)/values1(1,1);
    values13= eo_getconductorproperties('v3');
    C13= values13(1,2)/values1(1,1);
    values14 = eo_getconductorproperties('v4');
    C14= values14(1,2)/values1(1,1);

    ei_modifyconductorprop('v1',1,0);
    ei_modifyconductorprop('v2',1,1);
    ei_analyze(0)
    ei_loadsolution

    values2= eo_getconductorproperties('v2');
    C22= values2(1,2)/values2(1,1);
    values21 = eo_getconductorproperties('v1');
    C21= values21(1,2)/values2(1,1);
    values23= eo_getconductorproperties('v3');
    C23= values23(1,2)/values2(1,1);
    values24 = eo_getconductorproperties('v4');
    C24= values24(1,2)/values2(1,1);

    ei_modifyconductorprop('v2',1,0);
    ei_modifyconductorprop('v3',1,1);
    ei_analyze(0)
    ei_loadsolution

    values3= eo_getconductorproperties('v3');
    C33= values3(1,2)/values3(1,1);
    values31 = eo_getconductorproperties('v1');
    C31= values31(1,2)/values3(1,1);
    values32= eo_getconductorproperties('v2');
    C32= values32(1,2)/values3(1,1);
    values34 = eo_getconductorproperties('v4');
    C34= values34(1,2)/values3(1,1);

    ei_modifyconductorprop('v3',1,0);
    ei_modifyconductorprop('v4',1,1);
    ei_analyze(0)
    ei_loadsolution

    values4= eo_getconductorproperties('v4');
    C44= values4(1,2)/values4(1,1);
    values41 = eo_getconductorproperties('v1');
    C41= values41(1,2)/values4(1,1);
    values42= eo_getconductorproperties('v2');
    C42= values42(1,2)/values4(1,1);
    values43 = eo_getconductorproperties('v3');
    C43= values43(1,2)/values4(1,1);

    C = [C11,C12,C13,C14;C21,C22,C23,C24;C31,C32,C33,C34;C41,C42,C43,C44];


    %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       SIMULATION for L and Rac matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    openfemm
    newdocument(0);

    mi_getmaterial('Air');
    mi_getmaterial('18 AWG');
    mi_addmaterial('sheath_material', 1, 1, 0, 0, 35.33, 0, 0, 1, 0, 0, 0, 1, 1.02396529684335);
    mi_addmaterial('core', 1, 1, 0, 0, 58, 0, 0, 1, 0, 0, 0, 1, 1.02396529684335);
    mi_addmaterial('soil', 1, 1, 0, 0, 1e-8, 0, 0, 1, 0, 0, 0, 1, 0);

    %cable geometry-1
    %core
    mi_addnode(0.032,0);
    mi_addnode(-0.032,0);
    mi_addarc(0.032,0,-0.032,0,180,1);
    mi_addarc(-0.032,0,0.032,0,180,1);
    mi_addblocklabel(0,0); 
    mi_selectlabel(0,0);
    mi_setblockprop('core',0,0.001,'circuit1',0,0,1);
    mi_clearselected;

    %insulalation
    mi_addnode(0.0569,0);
    mi_addnode(-0.0569,0);
    mi_addarc(0.0569,0,-0.0569,0,180,1);
    mi_addarc(-0.0569,0,0.0569,0,180,1);
    mi_addblocklabel(-0.04445,0);
    mi_selectlabel(-0.04445,0);
    mi_setblockprop('Air',1,0,'Air',0,0,0);
    mi_clearselected;

    %sheath
    mi_addnode(0.0582,0);
    mi_addnode(-0.0582,0);
    mi_addarc(0.0582,0,-0.0582,0,180,1);
    mi_addarc(-0.0582,0,0.0582,0,180,1);
    mi_addblocklabel(0.05755,0); 
    mi_selectlabel(0.05755,0);
    mi_setblockprop('sheath_material',0,0.0005,'circuit2',0,0,1);
    mi_clearselected;

    % outer insulalation
    mi_addnode(0.0639,0);
    mi_addnode(-0.0639,0);
    mi_addarc(0.0639,0,-0.0639,0,180,1);
    mi_addarc(-0.0639,0,0.0639,0,180,1);
    mi_addblocklabel(-0.06105,0);
    mi_selectlabel(-0.06105,0);
    mi_setblockprop('Air',1,0,'Air',0,0,0);
    mi_clearselected;

    %cable geometry-2
    %core
    mi_addnode(0.468,0);
    mi_addnode(0.532,0);
    mi_addarc(0.468,0,0.532,0,180,1);
    mi_addarc(0.532,0,0.468,0,180,1);
    mi_addblocklabel(0.5,0); 
    mi_selectlabel(0.5,0);
    mi_setblockprop('core',0,0.001,'circuit3',0,0,1);
    mi_clearselected;

    %insulalation
    mi_addnode(0.4431,0);
    mi_addnode(0.5569,0);
    mi_addarc(0.4431,0,0.5569,0,180,1);
    mi_addarc(0.5569,0,0.4431,0,180,1);
    mi_addblocklabel(0.45555,0);
    mi_selectlabel(0.45555,0);
    mi_setblockprop('Air',1,0,'Air',0,0,0);
    mi_clearselected;

    %sheath
    mi_addnode(0.5582,0);
    mi_addnode(0.4418,0);
    mi_addarc(0.5582,0,0.4418,0,180,1);
    mi_addarc(0.4418,0,0.5582,0,180,1);
    mi_addblocklabel(0.55755,0); 
    mi_selectlabel(0.55755,0);
    mi_setblockprop('sheath_material',0,0.0005,'circuit4',0,0,1);
    mi_clearselected;

    % outer insulalation
    mi_addnode(0.5639,0);
    mi_addnode(0.4361,0);
    mi_addarc(0.5639,0,0.4361,0,180,1);
    mi_addarc(0.4361,0,0.5639,0,180,1);
    mi_addblocklabel(0.43895,0);
    mi_selectlabel(0.43895,0);
    mi_setblockprop('Air',1,0,'Air',0,0,0);
    mi_clearselected;

    %définition des conditions aux limites-boundary condition
    mi_addboundprop('A0',0,0,0,0,0,0,0,0,0);
    mi_drawline(-50,50,-50,-50);
    mi_drawline(50.5,50,50.5,-50);
    mi_drawline(-50,-50,50.5,-50);
    mi_drawline(-50,50,50.5,50);

    mi_selectsegment(-50,48);
    mi_setsegmentprop('A0',0,0,0,0);
    mi_clearselected;
    mi_selectsegment(-48,-50);
    mi_setsegmentprop('A0',0,0,0,0);
    mi_clearselected;
    mi_selectsegment(50.5,-48);
    mi_setsegmentprop('A0',0,0,0,0);
    mi_clearselected;
    mi_selectsegment(-48,50);
    mi_setsegmentprop('A0',0,0,0,0);
    mi_clearselected;

    mi_addblocklabel(0,5); %air
    mi_selectlabel(0,5);
    mi_setblockprop('Air',1,0,'Air',0,0,0);
    mi_clearselected;

    %soil line
    mi_drawline(-50,1.5,50.5,1.5);
    mi_addblocklabel(0.25,1);
    mi_selectlabel(0.25,1);
    mi_setblockprop('soil',1,0,'soil',0,0,0);
    mi_clearselected;

    mi_zoomnatural
    mi_saveas('.\femm\cableZ.FEM')

    %définition des circuits
    mi_addcircprop('circuit1',1,1); %core cable 1
    mi_addcircprop('circuit2',0,1); %sheath cable 1
    mi_addcircprop('circuit3',0,1); %core cable 2
    mi_addcircprop('circuit4',0,1); %sheath cable 2

    %i=1;
    for i = 1:max(size(freq))
       f = freq(1,i);

     mi_probdef(f,'meters','planar',1.e-8,l,30);
     mi_analyze
     mi_loadsolution

    % 1st row
        values1 = mo_getcircuitproperties('circuit1');
        L11(1,i)= real(values1(1,3))/values1(1,1);
        R11ac(1,i)= (1i)^2 *2*pi*f*(imag(values1(1,3))/values1(1,1));

        values12 = mo_getcircuitproperties('circuit2');
        L12(1,i)= real(values12(1,3))/values1(1,1);
        R12ac(1,i)=(1i)^2 *2*pi*f*(imag(values12(1,3))/values1(1,1));

        values13 = mo_getcircuitproperties('circuit3');
        L13(1,i)= real(values13(1,3))/values1(1,1);
        R13ac(1,i)= (1i)^2 *2*pi*f*(imag(values13(1,3))/values1(1,1));

        values14 = mo_getcircuitproperties('circuit4');
        L14(1,i)= real(values14(1,3))/values1(1,1);
        R14ac(1,i)=(1i)^2 *2*pi*f*(imag(values14(1,3))/values1(1,1));

        %modification des circuits
        mi_modifycircprop('circuit1',1,0);
        mi_modifycircprop('circuit2',1,1);

        mi_analyze
        mi_loadsolution

      %  2nd row
      values2 = mo_getcircuitproperties('circuit2');
        L22(1,i)= real(values2(1,3))/values2(1,1);
        R22ac(1,i)= (1i)^2 *2*pi*f*(imag(values2(1,3))/values2(1,1));

        values21 = mo_getcircuitproperties('circuit1');
        L21(1,i)= real(values21(1,3))/values2(1,1);
        R21ac(1,i)=(1i)^2 *2*pi*f*(imag(values21(1,3))/values2(1,1));

        values23 = mo_getcircuitproperties('circuit3');
        L23(1,i)= real(values23(1,3))/values2(1,1);
        R23ac(1,i)= (1i)^2 *2*pi*f*(imag(values23(1,3))/values2(1,1));

        values24 = mo_getcircuitproperties('circuit4');
        L24(1,i)= real(values24(1,3))/values2(1,1);
        R24ac(1,i)=(1i)^2 *2*pi*f*(imag(values24(1,3))/values2(1,1));

        %modification des circuits
        mi_modifycircprop('circuit2',1,0);
        mi_modifycircprop('circuit3',1,1);

        mi_analyze
        mi_loadsolution

     %  3rd row
        values3 = mo_getcircuitproperties('circuit3');
        L33(1,i)= real(values3(1,3))/values3(1,1);
        R33ac(1,i)= (1i)^2 *2*pi*f*(imag(values3(1,3))/values3(1,1));

        values31 = mo_getcircuitproperties('circuit1');
        L31(1,i)= real(values31(1,3))/values3(1,1);
        R31ac(1,i)=(1i)^2 *2*pi*f*(imag(values31(1,3))/values3(1,1));

        values32 = mo_getcircuitproperties('circuit2');
        L32(1,i)= real(values32(1,3))/values3(1,1);
        R32ac(1,i)= (1i)^2 *2*pi*f*(imag(values32(1,3))/values3(1,1));

        values34 = mo_getcircuitproperties('circuit4');
        L34(1,i)= real(values34(1,3))/values3(1,1);
        R34ac(1,i)=(1i)^2 *2*pi*f*(imag(values34(1,3))/values3(1,1));

        %modification des circuits
        mi_modifycircprop('circuit3',1,0);
        mi_modifycircprop('circuit4',1,1);

        mi_analyze
        mi_loadsolution

     %  4th row
      values4 = mo_getcircuitproperties('circuit4');
        L44(1,i)= real(values4(1,3))/values4(1,1);
        R44ac(1,i)= (1i)^2 *2*pi*f*(imag(values4(1,3))/values4(1,1));

        values41 = mo_getcircuitproperties('circuit1');
        L41(1,i)= real(values41(1,3))/values4(1,1);
        R41ac(1,i)=(1i)^2 *2*pi*f*(imag(values41(1,3))/values4(1,1));

        values43 = mo_getcircuitproperties('circuit3');
        L43(1,i)= real(values43(1,3))/values4(1,1);
        R43ac(1,i)= (1i)^2 *2*pi*f*(imag(values43(1,3))/values4(1,1));

        values42 = mo_getcircuitproperties('circuit2');
        L42(1,i)= real(values42(1,3))/values4(1,1);
        R42ac(1,i)=(1i)^2 *2*pi*f*(imag(values42(1,3))/values4(1,1));

        %modification des circuits
        mi_modifycircprop('circuit4',1,0);
        mi_modifycircprop('circuit1',1,1);

      Lf(:,:,i)= [L11(1,i),L12(1,i),L13(1,i),L14(1,i);L21(1,i),L22(1,i),L23(1,i),L24(1,i);L31(1,i),L32(1,i),L33(1,i),L34(1,i);L41(1,i),L42(1,i),L43(1,i),L44(1,i)];

      Rf(:,:,i)= Rdc+ [R11ac(1,i),R12ac(1,i),R13ac(1,i),R14ac(1,i);R21ac(1,i),R22ac(1,i),R23ac(1,i),R24ac(1,i);R31ac(1,i),R32ac(1,i),R33ac(1,i),R34ac(1,i);R41ac(1,i),R42ac(1,i),R43ac(1,i),R44ac(1,i)];

      %i=i+1;
    end

end