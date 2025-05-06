function [R_3D, R_XY, R_XZ, R_YZ] = Roundness_SVD(xyz_data,t_Qonset,t_Qoffset, centering) 
    %R_3D --> also referred to as "complexity"
    %R_plane where plane is subset of {XY,XZ,YZ} for loop projection

    Loop_data = xyz_data(:,t_Qonset:t_Qoffset);

    if strcmp(centering,'mean_centering')
        mean_coords = mean(Loop_data,2);
        Loop_data = Loop_data - mean_coords;
    elseif strcmp(centering,'mean_std_centering')
        Loop_data = normalize(Loop_data.', 'center', 'mean', 'scale');
        Loop_data = Loop_data.';
    end
    eig_values_3D = svd(Loop_data);
    R_3D = eig_values_3D(2) / eig_values_3D(1);
    
    %Projections
    x_loop = Loop_data(1,:);
    y_loop = Loop_data(2,:);
    z_loop = Loop_data(3,:);
    
    Projection_XY_data = [x_loop;y_loop];
    Projection_XZ_data = [x_loop;z_loop];
    Projection_YZ_data = [y_loop;z_loop];

    P_eigvalues = svd(Projection_XY_data);
    R_XY = P_eigvalues(2) / P_eigvalues(1);


    P_eigvalues = svd(Projection_XZ_data);

    R_XZ = P_eigvalues(2) / P_eigvalues(1);

    P_eigvalues = svd(Projection_YZ_data);

    R_YZ = P_eigvalues(2) / P_eigvalues(1);

end