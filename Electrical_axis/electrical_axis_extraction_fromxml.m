%% Setting Up
close all
clearvars
addpath /home/mka23/Desktop/UKBB/UKBB_Data/UKBB_ECGs_54K/
addpath dTools/
addpath fxns/
MyFolderInfo = dir('/home/mka23/Desktop/UKBB/UKBB_Data/UKBB_ECGs_54K/'); %%folder with all XML files (bulk data)

%% P1 - MAIN LOOP
ECG_Data = table(); % Initialize the table
for m = 3:length(MyFolderInfo)
    A=xml_read(MyFolderInfo(m).name);
    x88878_ID = str2double(MyFolderInfo(m).name(1:7));

    if isfield(A, 'RestingECGMeasurements') && isfield(A.RestingECGMeasurements, 'MedianSamples')
        [MedianBeat, s_Qonset, s_Qoffset] = MedianBeat_extract(A,false);
    else
        continue
    end
    xyz_data = leadcalc(MedianBeat([7:12,1,2],:),'kors'); %Using V1-6 + I & II (based on coefficient order in leadcalc(,kors)
    try
        [axis_vector_M1] = VCG_axis_maxNorm(xyz_data, s_Qonset, s_Qoffset,false);
        [axis_vector_M2] = VCG_axis_max_XYZ(xyz_data, s_Qonset, s_Qoffset,false);
        [axis_vector_M3] = VCG_axis_average(xyz_data, s_Qonset, s_Qoffset,false);
        [axis_vector_M4] = VCG_axis_VW_average(xyz_data, s_Qonset, s_Qoffset,false);
        [axis_vector_M5, ~, ~] = VCG_axis_max_SVD(xyz_data, s_Qonset, s_Qoffset, false, 'mean_centering'); %can change centering here
    catch
        continue
    end
    %Unit vectors
    axis_vector_M1 = axis_vector_M1/norm(axis_vector_M1);
    axis_vector_M2 = axis_vector_M2/norm(axis_vector_M2);
    axis_vector_M3 = axis_vector_M3/norm(axis_vector_M3);
    axis_vector_M4 = axis_vector_M4/norm(axis_vector_M4);
    axis_vector_M5 = axis_vector_M5/norm(axis_vector_M5);
    
    RR = A.RestingECGMeasurements.RRInterval.CONTENT;
    Pd = A.RestingECGMeasurements.PDuration.CONTENT;
    PQ = A.RestingECGMeasurements.PQInterval.CONTENT; 
    QRSd = A.RestingECGMeasurements.QRSDuration.CONTENT; 
    QT = A.RestingECGMeasurements.QTInterval.CONTENT;
    QTc = A.RestingECGMeasurements.QTCInterval.CONTENT;
    Paxis = A.RestingECGMeasurements.PAxis.CONTENT;
    Raxis = A.RestingECGMeasurements.RAxis.CONTENT;
    Taxis = A.RestingECGMeasurements.TAxis.CONTENT;
    

    duration_values = {RR, Pd, PQ, QRSd, QT, QTc, Paxis, Raxis, Taxis};
    for i = 1:numel(duration_values)
        if isempty(duration_values{i})
            duration_values{i} = NaN;
        end
    end
    % Reassign the values to their original variables
    RR = duration_values{1}; Pd = duration_values{2}; PQ = duration_values{3};
    QRSd = duration_values{4}; QT = duration_values{5}; QTc = duration_values{6};
    Paxis = duration_values{7}; Raxis = duration_values{8}; Taxis = duration_values{9};

    t_Qonset = A.RestingECGMeasurements.QOnset.CONTENT;
    t_Qoffset = A.RestingECGMeasurements.QOffset.CONTENT;
    t_Toffset = A.RestingECGMeasurements.TOffset.CONTENT; 
    t_Ponset = A.RestingECGMeasurements.POnset.CONTENT;
    t_Poffset = A.RestingECGMeasurements.POffset.CONTENT;

    % converting [] to NaN
    time_values = {t_Ponset, t_Poffset, t_Qonset, t_Qoffset, t_Toffset};
    
    for i = 1:numel(time_values)
        if isempty(time_values{i})
            time_values{i} = NaN;
        end
    end
    
    % Reassign the values 
    t_Ponset = time_values{1}; t_Poffset = time_values{2};
    t_Qonset = time_values{3}; t_Qoffset = time_values{4}; t_Toffset = time_values{5};

   
    QRS_area_1 = Calc_area(xyz_data,s_Qonset,s_Qoffset,'mean_O');
    Roundness_1 = QRS_area_1/norm(axis_vector_M1)^2;
    

    [R_3D, R_XY, R_XZ, R_YZ] = Roundness_SVD(xyz_data,s_Qonset,s_Qoffset, 'mean_centering'); 
    [Pl_alpha, Pl_beta, Pl_gamma] = planarity_SVD(xyz_data,s_Qonset,s_Qoffset, 'mean_centering');
    [RMSE, R_squared,RMSE_perp] = LSE_planarity(xyz_data, s_Qonset,s_Qoffset, false);
    Pl_R2 = R_squared;
    
    dipole_mag = [];
    for k = 1: length(xyz_data)
        dipole_mag(k) = norm(xyz_data(:,k));
    end

    if isfield(A, 'Interpretation') && isfield(A.Interpretation, 'Conclusion') && ...
            isfield(A.Interpretation.Conclusion, 'ConclusionText') && ~isempty(A.Interpretation.Conclusion.ConclusionText)
        % Extract the ConclusionText if it exists and is not empty
        conclusion = strjoin(A.Interpretation.Conclusion.ConclusionText, ', ');
    else
        % Handle the case where ConclusionText is empty or missing
        conclusion = 'Conclusion not found or is empty';
    end    
    
    

    %%EVERYTHING:
    patientRow = table(x88878_ID,...
        {mat2str(axis_vector_M1)},{mat2str(axis_vector_M2)},{mat2str(axis_vector_M3)},{mat2str(axis_vector_M4)},{mat2str(axis_vector_M5)},... 
        {MedianBeat},{xyz_data}, PQ, Pd, QRSd, QT, QTc, RR, t_Ponset,t_Poffset, t_Qonset, t_Qoffset, t_Toffset,...
        Paxis, Raxis, Taxis, {dipole_mag},...
        QRS_area_1, Roundness_1, R_3D, R_XY, R_XZ, R_YZ, ...
        Pl_R2, Pl_alpha, Pl_beta, Pl_gamma,{conclusion},...
        'VariableNames', ...
        {'x88878_ID', ...
        'axis_vector_M1','axis_vector_M2','axis_vector_M3','axis_vector_M4','axis_vector_M5',... 
        'MedianBeat','xyz_data', 'PQ', 'Pd', 'QRSd', 'QT', 'QTc', 'RR', 't_Ponset', 't_Poffset', 't_Qonset', 't_Qffset', 't_Toffset',...
        'Paxis', 'Raxis', 'Taxis', 'dipole_mag',...
        'QRS_area_1', 'Roundness_1',  'R_3D', 'R_XY', 'R_XZ', 'R_YZ', ...
        'Pl_R2', 'Pl_alpha', 'Pl_beta', 'Pl_gamma', 'conclusion'});

    disp(m)
    % Append each patient data row to the table
    ECG_Data = [ECG_Data; patientRow];
end


VCG_vectors = ECG_Data(:,1:6);
writetable(VCG_vectors, ['Extracted_raw_data/',num2str(height(ECG_Data)), ...
   '_VCG_vectors.csv']);

ECG_metrics = ECG_Data(:,[1,9:22,34]);
writetable(ECG_metrics, ['Extracted_raw_data/',num2str(height(ECG_Data)), ...
   '_ECG_metrics.csv'],'Delimiter',',','QuoteStrings','all');

VCG_metrics = ECG_Data(:,[1,24:33]);
writetable(VCG_metrics, ['Extracted_raw_data/',num2str(height(ECG_Data)), ...
   '_VCG_metrics.csv'],'Delimiter',',','QuoteStrings','all');

time_series_MedianBeat = ECG_Data(:,7);
writetable(time_series_MedianBeat, ['Extracted_raw_data/',num2str(height(ECG_Data)), ...
   '_ECG_12lead.csv'],'Delimiter',',','QuoteStrings','all');



