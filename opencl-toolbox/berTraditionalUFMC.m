%% Configure Repeatability State
% s = rng(211);                                           % Set RNG state for repeatability

%% Initial Parameters
numFFT = 512;                                           % number of FFT points
subbandSize = 20;                                       % must be > 1 
numSubbands = 10;                                       % numSubbands*subbandSize <= numFFT
subbandOffset = numFFT/2-subbandSize*numSubbands/2;     %band center in ofdm block
filterLen = 43;                                         % similar to cyclic prefix length
slobeAtten = 40;                                        % side-lobe attenuation, dB
bitsPerSubCarrier = 4;                                  % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
snrdB = 15;                                             % SNR in dB
prototypeFilter = chebwin(filterLen, slobeAtten);       % Design window with specified attenuation
errorRate = comm.ErrorRate("ResetInputPort",true);      % BER System Object
EbNoVec = (0:15)';
snrVec = EbNoVec + 10*log10(bitsPerSubCarrier) + 10*log10(numSubbands*subbandSize/numFFT);
berVec = zeros(length(EbNoVec),3);
errorStats = zeros(1,3);
maxBitErrors = 400;    % Maximum number of bit errors
maxNumBits = 1e7;      % Maximum number of bits transmitted

for m = 1:length(EbNoVec)
    snr = snrVec(m);

    while errorStats(2) <= maxBitErrors && errorStats(3) <= maxNumBits
    
        %% Transmit-end processing
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
        
        %% Receive-end Processing
        % channel = comm.AWGNChannel('NoiseMethod','Variance', 'VarianceSource','Input port');    % Noise Channel
        % powerDB = 10*log10(var(txSig));                                                         % Calculate Tx signal power
        % noiseVar = 10.^(0.1*(powerDB-snrdB));                                                   % Calculate the noise variance
        % rxSig = channel(txSig,noiseVar);                                                        % Pass the signal through a noisy channel
        rxSig = awgn(txSig, snr, 'measured');
        yRxPadded = [rxSig; zeros(2*numFFT-numel(txSig),1)];                                    % zero padding to 2*FFT lenght
        RxSymbols2x = fftshift(fft(yRxPadded));                                                 % perform FFT
        RxSymbols = RxSymbols2x(1:2:end);                                                       % downsample by 2
        dataRxSymbols = RxSymbols(subbandOffset+(1:numSubbands*subbandSize));                   % Select data subcarriers
        
        % Use zero-forcing equalizer after OFDM demodulation
        rxf = [prototypeFilter.*exp(1i*2*pi*0.5*(0:filterLen-1)'/numFFT); ...
               zeros(numFFT-filterLen,1)];
        prototypeFilterFreq = fftshift(fft(rxf));
        prototypeFilterInv = 1./prototypeFilterFreq(numFFT/2-subbandSize/2+(1:subbandSize));

        % Equalize per subband - undo the filter distortion
        dataRxSymbolsMat = reshape(dataRxSymbols,subbandSize,numSubbands);
        EqualizedRxSymbolsMat = bsxfun(@times,dataRxSymbolsMat,prototypeFilterInv);
        EqualizedRxSymbols = EqualizedRxSymbolsMat(:);
        
        rxBits = qamdemod(EqualizedRxSymbols, 2^bitsPerSubCarrier, 'OutputType', 'bit', ...
            'UnitAveragePower', true);
        
        errorStats = errorRate(inpData(:), rxBits, 0);
    end

    berVec(m,:) = errorStats;                              % Save BER data
    errorStats = errorRate(inpData(:), rxBits, 1);         % Reset the error rate calculator
end

%%
berTheory = berawgn(EbNoVec,'qam',2^bitsPerSubCarrier,'nondiff');
% figure (1)
semilogy(EbNoVec,berVec(:,1), '--o')
% hold on
% semilogy(EbNoVec,berTheory)
% legend('Simulation','Theory','Location','Best')
% xlabel('Eb/No (dB)')
% ylabel('Bit Error Rate')
% grid on
% hold off