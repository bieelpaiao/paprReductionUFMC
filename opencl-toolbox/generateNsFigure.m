clear; close all; clc;

numFFT_values = [64 128 256 512 1024 2048]; % vector of number of FFT points
numSubbands_values = [2 4 10 20 40 80]; % vector of numSubbands (numSubbands*subbandSize <= numFFT)

% Dolph-Chebyshev window design parameters
filterLen = 43;      % similar to cyclic prefix length
slobeAtten = 40;     % side-lobe attenuation, dB

bitsPerSubCarrier = 4;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
snrdB = 15;              % SNR in dB

% Design window with specified attenuation
prototypeFilter = chebwin(filterLen, slobeAtten);

beta = 0.1;

% Inicialização do gráfico CCDF
figure;
colors = lines(length(numFFT_values)); % Para usar cores diferentes para cada curva

for i = 1:length(numFFT_values)
    rng(211);       % Set RNG state for repeatability
    
    numFFT = numFFT_values(i); % number of FFT points
    numSubbands = numSubbands_values(i); % must be > 1 
    subbandSize = 20;    % must be > 1 

    % Transmit-end processing
    %  Initialize arrays
    inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands);
    txSig = complex(zeros(numFFT+filterLen-1, 1));
    
    nOFDM = numSubbands*subbandSize;
    
    np = round(nOFDM*beta);
    
    pim = matrixGeneration(nOFDM, np);
    
    bits = randi([0 1], bitsPerSubCarrier*subbandSize*numSubbands, 1);
    symbols = qammod(bits, 2^bitsPerSubCarrier, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);
    
    symbolsPrecoded = precod(pim, symbols);
    
    subbandSize = length(symbolsPrecoded)/numSubbands;
    
    subbandOffset = numFFT/2-subbandSize*numSubbands/2;
    
    S2P = reshape(symbolsPrecoded, subbandSize, numSubbands);
    
    % Loop over each subband
    for bandIdx = 1:numSubbands
    
        symbolsIn = S2P(:, bandIdx);
    
        % Pack subband data into an OFDM symbol
        offset = subbandOffset+(bandIdx-1)*subbandSize; 
        symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
                         zeros(numFFT-offset-subbandSize, 1)];
        ifftOut = ifft(ifftshift(symbolsInOFDM));
    
        % Filter for each subband is shifted in frequency
        bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT* ...
                     ((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );    
        filterOut = conv(bandFilter,ifftOut);
    
        % Sum the filtered subband responses to form the aggregate transmit
        % signal
        txSig = txSig + filterOut;     
    end
    
    % Compute peak-to-average-power ratio (PAPR)
    pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);
    paprUFMC = pm(txSig);
    disp(['Peak-to-Average-Power-Ratio (PAPR) for UFMC = ' num2str(paprUFMC) ' dB']);
    plotCCDF(pm, 'LineWidth', 2);
    if (i == 1)
        hold on
    end
end

hold off;
title('')
xlabel('$PAPR_{0}$', 'interpreter', 'latex')
ylabel('$P(PAPR \geq PAPR_{0})$','interpreter','latex')
legend('N=64, n=2','N=128, n=4','N=256, n=10','N=512, n=20','N=1024, n=40','N=2048, n=80')
set(gca, 'LineWidth', 1.5);
grid on

%% Save Image
% Load Style Sheet created on Figure > File > Export Setup
% OBS: hgexport will not be supported in a future release
s = hgexport('readstyle','IC');
% s.Format = 'png';
s.Format = 'eps';

% Export image as .png and .fig
% hgexport(gcf,'./img/png/UFMC_CCDF_Ns.png',s);
hgexport(gcf,'./img/eps/UFMC_CCDF_Ns.eps',s);
% savefig("./img/fig/UFMC_CCDF_Ns.fig")

