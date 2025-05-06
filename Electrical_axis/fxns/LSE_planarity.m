function [RMSE, R_squared,RMSE_perp] = LSE_planarity(xyz_data, start_time,end_time, PLOT);

    %Plane in the form: z= ax + by + c --> Matrix form: A*fit = B 
    x = xyz_data(1, start_time:end_time);
    y = xyz_data(2, start_time:end_time); 
    z = xyz_data(3, start_time:end_time);
    
    % Fitting
    tmp_A = [];
    tmp_b = [];
    
    for i = 1:length(x)
        tmp_A = [tmp_A; x(i), y(i), 1];
        tmp_b = [tmp_b; z(i)];
    end
    
    A = tmp_A;
    B = tmp_b;
    
    % Manual solution
    fit = (A' * A) \ (A' * B);
    a=fit(1); b=fit(2);c=fit(3);
    errors = B - A * fit;
    
    Z_plane = A * fit;   % Calculate Surface At Each Data Triplet
    
    % Plot the plane, and points with error visible
    if PLOT == true
        [X, Y] = meshgrid(linspace(min(x), max(x), 10), linspace(min(y), max(y), 10));
        Z = a * X + b * Y + c;
        figure;
        surf(X, Y, Z, 'FaceAlpha', 0.5);
        hold on;
        scatter3(x, y, z, 'filled', 'b');
        hold on;
        plot3([x; x], [y; y], [z; Z_plane'], '-r', 'LineWidth',1)  % Plot Errors (Red Lines From Plane To Data)
        xlabel('x');
        ylabel('y');
        zlabel('z');
        title('Best-Fit Plane for QRS Loop Data points');
    end
    
    % sum of squared residuals (SSR)
    ssr = sum(errors.^2);
    % total sum of squares (SST)
    mean_z = mean(z);
    sst = sum((z - mean_z).^2);
    
    % Calculate R-squared (R²)
    R_squared = 1 - ssr / sst;
    RMSE = sqrt(mean(errors.^2));
    
    
    perp_dist = [];
    for i = 1:(end_time-start_time)
    perp_dist = [perp_dist; abs(a* x(i) + b * y(i) - z(i) + c ) / sqrt(fit(1)^2 + fit(2)^2 + 1^2)];
    end
    
    RMSE_perp = sqrt(mean(perp_dist.^2));

end