clear; close all; clc;

s = rng(211);

beta_values = [0 0.1 0.25 0.5 0.75]; % vector of beta values

% Inicialização do gráfico CCDF
figure;
colors = lines(length(beta_values)); % Para usar cores diferentes para cada curva
  
for i = 1:length(beta_values)
    s = rng(211);       % Set RNG state for repeatability

    beta = beta_values(i);

    if beta == 0 % without precoding
        numFFT = 512;                                           % number of FFT points
        subbandSize = 20;                                       % must be > 1 
        numSubbands = 10;                                       % numSubbands*subbandSize <= numFFT
        subbandOffset = numFFT/2-subbandSize*numSubbands/2;     %band center in ofdm block
        filterLen = 43;                                         % similar to cyclic prefix length
        slobeAtten = 40;                                        % side-lobe attenuation, dB
        bitsPerSubCarrier = 4;                                  % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
        snrdB = 15;                                             % SNR in dB
        prototypeFilter = chebwin(filterLen, slobeAtten);       % Design window with specified attenuation

        %  Initialize arrays
        inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands);
        txSig = complex(zeros(numFFT+filterLen-1, 1));

        %  Loop over each subband
        for bandIdx = 1:numSubbands
        
            bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize, 1);
            % QAM Symbol mapper
            symbolsIn = qammod(bitsIn, 2^bitsPerSubCarrier, 'InputType', 'bit', ...
            'UnitAveragePower', true);
            inpData(:,bandIdx) = bitsIn; % log bits for comparison
            
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
        hold on
    else
        numFFT = 512; % number of FFT points
        numSubbands = 10; % numSubbands*subbandSize <= numFFT
        subbandSize = 20;    % must be > 1
        
        % Dolph-Chebyshev window design parameters
        filterLen = 43;      % similar to cyclic prefix length
        slobeAtten = 40;     % side-lobe attenuation, dB
        
        bitsPerSubCarrier = 4;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
        snrdB = 15;              % SNR in dB
        
        % Design window with specified attenuation
        prototypeFilter = chebwin(filterLen, slobeAtten);

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
    end
end

hold off;
title('')
xlabel('$PAPR_{0}$', 'interpreter', 'latex')
ylabel('$P(PAPR \geq PAPR_{0})$','interpreter','latex')
legend('Traditional','Precoding (β=10%)','Precoding (β=25%)','Precoding (β=50%)','Precoding (β=75%)')
set(gca, 'LineWidth', 1.5);
grid on

%% Save Image
% Load Style Sheet created on Figure > File > Export Setup
% OBS: hgexport will not be supported in a future release
s = hgexport('readstyle','IC');
% s.Format = 'png';
s.Format = 'eps';

% Export image as .png and .fig
% hgexport(gcf,'./img/png/UFMC_CCDF_Betas.png',s);
hgexport(gcf,'./img/eps/UFMC_CCDF_Betas.eps',s);
% savefig("./img/fig/UFMC_CCDF_Betas.fig")
