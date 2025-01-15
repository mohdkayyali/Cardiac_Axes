function [MedianBeat, t_Qonset, t_Qoffset] = MedianBeat_extract(ECG_DATA_XML,PLOT)
    A = ECG_DATA_XML; %Struct file
    
    % Median Beat reading
    MedianBeat= [];
    Res = A.RestingECGMeasurements.MedianSamples.Resolution.CONTENT;
    Fs=A.RestingECGMeasurements.MedianSamples.SampleRate.CONTENT;
    t_end = A.RestingECGMeasurements.MedianSamples.LastValid.CONTENT;
    
    %All leads
    for k=1:12
        lead_points = (Res/1000) * cell2mat(textscan(A.RestingECGMeasurements.MedianSamples.WaveformData(k).CONTENT,'%f','Delimiter',','))';
        MedianBeat=[ MedianBeat; lead_points(1:t_end)];
    end

    if PLOT==true
         figure(),plot((0: 1:t_end-1) / Fs,MedianBeat),hold on
         title('Median Beat:  all leads'), grid on
         legend('I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'), xlabel('Time (seconds)'), ylabel('Amplitude (mV)');
    end

    t_Qonset = A.RestingECGMeasurements.QOnset.CONTENT * 0.5;
    t_Qoffset = A.RestingECGMeasurements.QOffset.CONTENT * 0.5;
   
 end


