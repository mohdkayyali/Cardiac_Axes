function [axis_vector] = VCG_axis_VW_average(xyz_data, t_Qonset, t_Qoffset, PLOT)
    
    Loop_data = xyz_data(:,t_Qonset:t_Qoffset);
    Fs = 500; %Sampling frequency
    T = 1/Fs;
    
    dP = diff(Loop_data, 1, 2);
    dist_3d = sqrt(sum(dP.^2, 1)); 
    
    % Calculate the velocity as the rate of change of 3D distance
    velocity_data = dist_3d / T;
        
    weightings = velocity_data./max(velocity_data);
    VWA = sum(Loop_data(:,2:length(Loop_data)) .* weightings,2);

    axis_vector = VWA;

    QRS_start = xyz_data(:,t_Qonset + 1); % +1 for true isolelectric start

    % 
    if PLOT == true 
    figure(), plot3(Loop_data(1,:),Loop_data(2,:),Loop_data(3,:)), xlabel('X Amplitude (mV)'), ylabel('Y Amplitude (mV)'), zlabel('Z Amplitude (mV)'), hold on
    scatter3(QRS_start(1), QRS_start(2), QRS_start(3), 'red', 'filled');  hold on
    scatter3(axis_vector(1), axis_vector(2),axis_vector(3), 'green', 'filled'), hold on
    plot3([QRS_start(1) axis_vector(1)],[QRS_start(2) axis_vector(2)],[QRS_start(3) axis_vector(3)])
    legend('QRS Loop','QRS start','Point of mean XYZ','Electrical Axis')
    grid on
    rotate3d on

    % figure(), subplot(311),plot(xyz_data(1,:)), hold on, yline(axis_vector(1),'red'), ylabel('X Amplitude (mV'), grid on,...
    %             subplot(312),plot(xyz_data(2,:)), hold on, yline(axis_vector(2)),ylabel('Y Amplitude (mV'),grid on, ...
    %             subplot(313),plot(xyz_data(3,:)), hold on, yline(axis_vector(3)),ylabel('Z Amplitude (mV'),grid on;
end