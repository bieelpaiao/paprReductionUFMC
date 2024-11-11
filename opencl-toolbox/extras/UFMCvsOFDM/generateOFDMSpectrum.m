clear; close all; clc;

s = rng(211);       % Set RNG state for repeatability

numFFT = 512;        % number of FFT points
subbandSize = 20;    % must be > 1 
numSubbands = 10;    % numSubbands*subbandSize <= numFFT
subbandOffset = 156; % numFFT/2-subbandSize*numSubbands/2 for band center

bitsPerSubCarrier = 4;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
snrdB = 15;              % SNR in dB

% Transmit-end processing
%  Initialize arrays
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands);

for bandIdx = 1:numSubbands
    bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize, 1);
    inpData(:,bandIdx) = bitsIn;
end

symbolsIn = qammod(inpData(:), 2^bitsPerSubCarrier, 'InputType', 'bit', ...
    'UnitAveragePower', true);

% Process all sub-bands together
offset = subbandOffset; 
symbolsInOFDM = [zeros(offset, 1); symbolsIn; ...
                 zeros(numFFT-offset-subbandSize*numSubbands, 1)];
ifftOut = sqrt(numFFT).*ifft(ifftshift(symbolsInOFDM));

% Plot power spectral density (PSD) over all subcarriers
[psd,f] = periodogram(ifftOut, rectwin(length(ifftOut)), numFFT*2, ...
                      1, 'centered'); 
hFig1 = figure; 
plot(f,10*log10(psd)); 
grid on
axis([-0.5 0.5 -100 20]);
xlabel('Normalized frequency'); 
ylabel('PSD (dBW/Hz)')
title(['OFDM, ' num2str(numSubbands*subbandSize) ' Subcarriers'])
set(hFig1, 'Position', figposition([46 50 25 30]));
set(gca, 'LineWidth', 1.5);

%% Save Image
% Load Style Sheet created on Figure > File > Export Setup
% OBS: hgexport will not be supported in a future release
s = hgexport('readstyle','UFMCvsOFDM');
s.Format = 'png';

% Export image as .png and .fig
hgexport(gcf,'./img/png/OFDMSpectrum.png',s);
savefig('./img/fig/OFDMSpectrum.fig')
