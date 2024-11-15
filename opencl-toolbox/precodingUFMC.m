s = rng(211);       % Set RNG state for repeatability

numFFT = 512;        % number of FFT points
subbandSize = 20;    % must be > 1 
numSubbands = 10;    % numSubbands*subbandSize <= numFFT

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
beta = 0.35;
np = round(nOFDM*beta);

pim = matrixGeneration(nOFDM, np);

bits = randi([0 1], bitsPerSubCarrier*subbandSize*numSubbands, 1);
symbols = qammod(bits, 2^bitsPerSubCarrier, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);

symbolsPrecoded = precodCL(pim, symbols);

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
% plotCCDF(pm);
