function [axis_vector, eigenvalues, eigenvectors] = VCG_axis_max_SVD(xyz_data, t_Qonset, t_Qoffset, PLOT, centering)
    Loop_data = xyz_data(:,t_Qonset:t_Qoffset);

    if strcmp(centering,'mean_centering')
        mean_coords = mean(Loop_data,2);
        Loop_data = Loop_data - mean_coords;
    elseif strcmp(centering,'mean_std_centering')
        Loop_data = normalize(Loop_data.', 'center', 'mean', 'scale');
        Loop_data = Loop_data.';
    end
    %Singular Value Decomposition (SVD)
    [U, S, l] = svd(Loop_data);
    
    % Extract singular values "eigenvalues" with corresponding left singular
    % vectors "eigenvectors"
    eigenvalues = diag(S);
    eigenvectors = U;
    v1 = U(:, 1);
    % lambda2 = singular_values(2); v2 = U(:, 2);
    % lambda3 = singular_values(3); v3 = U(:, 3);
    
    axis_vector = v1;
    
    if PLOT == true

        x_QRS = Loop_data(1,:);
        y_QRS = Loop_data(2,:);
        z_QRS = Loop_data(3,:);
        QRS_start = xyz_data(:,t_Qonset+1);

        figure()
        plot3(x_QRS, y_QRS, z_QRS, linewidth = 1.5); hold on
        scatter3(x_QRS, y_QRS, z_QRS); hold on
        plot3([QRS_start(1) v1(1)],[QRS_start(2) v1(2)],[QRS_start(3) v1(3)], LineWidth=2); hold on
        plot3([QRS_start(1) axis_vector_M1(1)],[QRS_start(2) axis_vector_M1(2)],[QRS_start(3) axis_vector_M1(3)], LineWidth=2);
        xlabel('X Amplitude (mV)'); ylabel('Y Amplitude (mV)'); zlabel('Z Amplitude (mV)');
        title('3D QRS Loop'); legend('QRS Loop','QRS Vector: Max norm', 'QRS Vector: Max SV'); grid on; 
    end


end