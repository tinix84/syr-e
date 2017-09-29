function [dataSet] = SaveInformation(geo,mat,dataSet)
%% ========================================================================
%% This function is used to make order in data saved in temporary folders 
%%
%% ========================================================================

dataSet.AirGapThickness = geo.g; % airgap thickness
dataSet.AirGapThickness = round(dataSet.AirGapThickness,2);
dataSet.AirGapRadius =    geo.r; % machine airgap radius
dataSet.AirGapRadius =    round(dataSet.AirGapRadius,2);
dataSet.ToothLength =     geo.lt; % tooth length
dataSet.ToothLength =     round(dataSet.ToothLength,2);
dataSet.StatorSlotOpen =  geo.acs; % stator slot open in [p.u.]
dataSet.StatorSlotOpen =  round(dataSet.StatorSlotOpen,2);
dataSet.ToothWidth =      geo.wt; % Bgap/Btooth (determines tooth width^-1, yoke width^-1)
dataSet.ToothWidth =      round(dataSet.ToothWidth,2);
dataSet.ToothTangDepth =  geo.ttd; % tooth tang depth [mm]
dataSet.ToothTangDepth =  round(dataSet.ToothTangDepth,2);
dataSet.Br =              round(mat.LayerMag.Br,4); % Br
dataSet.Br =              round(dataSet.Br,2);
dataSet.ALPHApu =         geo.dalpha_pu;
dataSet.ALPHApu =         round(dataSet.ALPHApu,2);
dataSet.HCpu =            geo.hc_pu;
dataSet.HCpu =            round(dataSet.HCpu,2);
dataSet.DepthOfBarrier =  geo.dx;    % the depth of the barriers radial-wise in per unit or number of segments for SPM
dataSet.DepthOfBarrier =  round(dataSet.DepthOfBarrier,2);
dataSet.RQ =              geo.RQ;


end

