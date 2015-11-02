function [ho,mag,phase] = FFTAnalysis(y,maxHO,plt)

%*** FFTAnalysis (v1.00) **************************************************
% [ho,mag,thd] = FFTAnalysis(y,maxHO,plt)
% 
% INPUT:
%        y: segnale da analizzare
%    maxHO: massimo ordine armonico rappresentato nel grafico
%      plt: flag di abilitazione del plot
%
% OUTPUT:
%       ho: harmonic order
%      mag: modulo FFT
% *************************************************************************
% S. Ettorre
% *************************************************************************

maxIdx = maxHO+1;
ho = [1:maxIdx]-1;

% Performs the FFT analysis
nPts = length(y);
yFFT = fft(y,nPts);

% Calcola ampiezze
mag = 2*abs(yFFT)/nPts;
mag(1) = mag(1)/2;
mag = mag(1:maxIdx);

phase = angle(yFFT);
phase = phase(1:maxIdx);

if plt==1
    magPlot = mag;
    hoPlot = ho;
    bar(hoPlot,magPlot',0.1);
    grid
    xlabel('Harmonic Order');
    ylabel('Magnitude');    
end
