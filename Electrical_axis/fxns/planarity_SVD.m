function [alpha, beta, gamma] = planarity_SVD(xyz_data,t_Qonset,t_Qoffset, centering ) 
    %alpha, beta, gamma are 3 measure of planarity
    %alpha --> [0,1] 0=planar, 1=extends equally in 3 directions
    %beta --> 1 = Planar
    %gamma --> 0= Planar
    Loop_data = xyz_data(:,t_Qonset:t_Qoffset);

    if strcmp(centering,'mean_centering')
        mean_coords = mean(Loop_data,2);
        Loop_data = Loop_data - mean_coords;
    elseif strcmp(centering,'mean_std_centering')
        Loop_data = normalize(Loop_data.', 'center', 'mean', 'scale');
        Loop_data = Loop_data.';
    end
    % Extract singular values "eigenvalues" with corresponding left singular
    % vectors "eigenvectors"
    eigenvalues =svd(Loop_data);


    alpha = eigenvalues(3) / eigenvalues(1);
    
    beta = sum(eigenvalues(1:2)) / sum(eigenvalues);

    gamma = eigenvalues(3) / sum(eigenvalues(1:2));

end