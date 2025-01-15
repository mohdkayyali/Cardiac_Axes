function [axis_vector] = VCG_axis_max_XYZ(xyz_data, t_Qonset, t_Qoffset, PLOT)
    QRS_start = xyz_data(:,t_Qonset + 1); % +1 for true isolelectric start
    Loop_data = xyz_data(:,t_Qonset:t_Qoffset);

    %Find the maximum magnitudes of the components in raw values
    [max_magnitude, max_index] = max(abs(Loop_data), [], 2);
    x_index = max_index(1);
    y_index = max_index(2);
    z_index = max_index(3);

    % Create a new vector with components having the largest magnitudes in raw values
    axis_vector = [Loop_data(1, x_index); Loop_data(2, y_index);Loop_data(3, z_index)];
    % 
    if PLOT == true

    figure(), plot3(Loop_data(1,:),Loop_data(2,:),Loop_data(3,:)), xlabel('X Amplitude (mV)'), ylabel('Y Amplitude (mV)'), zlabel('Z Amplitude (mV)'), hold on
    scatter3(QRS_start(1), QRS_start(2), QRS_start(3), 'red', 'filled');  hold on
    scatter3(axis_vector(1), axis_vector(2),axis_vector(3), 'green', 'filled'), hold on
    plot3([QRS_start(1) axis_vector(1)],[QRS_start(2) axis_vector(2)],[QRS_start(3) axis_vector(3)])
    legend('QRS Loop','QRS start','Point of max XYZ','Electrical Axis')
    grid on
    rotate3d on

    figure(), subplot(311),plot(xyz_data(1,:)), hold on, scatter(x_index+t_Qonset-1,axis_vector(1),'red'), ylabel('X Amplitude (mV'), grid on,...
                subplot(312),plot(xyz_data(2,:)), hold on, scatter(y_index+t_Qonset-1,axis_vector(2)),ylabel('Y Amplitude (mV'),grid on, ...
                subplot(313),plot(xyz_data(3,:)), hold on, scatter(z_index+t_Qonset-1,axis_vector(3)),ylabel('Z Amplitude (mV'),grid on;
    end
end


