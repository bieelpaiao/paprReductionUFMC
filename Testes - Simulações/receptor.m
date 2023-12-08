% % Add WGN
% rxSig = awgn(txSig, snrdB, 'measured');
% 
% % Pad receive vector to twice the FFT Length (note use of txSig as input)
% %   No windowing or additional filtering adopted
% yRxPadded = [rxSig; zeros(2*numFFT-numel(txSig),1)];
% 
% % Perform FFT and downsample by 2
% RxSymbols2x = fftshift(fft(yRxPadded));
% RxSymbols = RxSymbols2x(1:2:end);
% 
% % Select data subcarriers
% dataRxSymbols = RxSymbols(subbandOffset+(1:numSubbands*subbandSize));
% 
% % Plot received symbols constellation
% % constDiagRx = comm.ConstellationDiagram('ShowReferenceConstellation', ...
% %     false, 'Position', figposition([20 15 25 30]), ...
% %     'Title', 'UFMC Pre-Equalization Symbols', ...
% %     'Name', 'UFMC Reception', ...
% %     'XLimits', [-150 150], 'YLimits', [-150 150]);
% % constDiagRx(dataRxSymbols);
% 
% % Use zero-forcing equalizer after OFDM demodulation
% rxf = [prototypeFilter.*exp(1i*2*pi*0.5*(0:filterLen-1)'/numFFT); ...
%        zeros(numFFT-filterLen,1)];
% prototypeFilterFreq = fftshift(fft(rxf));
% prototypeFilterInv = 1./prototypeFilterFreq(numFFT/2-subbandSize/2+(1:subbandSize));
% 
% % Equalize per subband - undo the filter distortion
% dataRxSymbolsMat = reshape(dataRxSymbols,subbandSize,numSubbands);
% EqualizedRxSymbolsMat = bsxfun(@times,dataRxSymbolsMat,prototypeFilterInv);
% EqualizedRxSymbols = EqualizedRxSymbolsMat(:);
% 
% % Plot equalized symbols constellation
% % constDiagEq = comm.ConstellationDiagram('ShowReferenceConstellation', ...
% %     false, 'Position', figposition([46 15 25 30]), ...
% %     'Title', 'UFMC Equalized Symbols', ...
% %     'Name', 'UFMC Equalization');
% % constDiagEq(EqualizedRxSymbols);
% 
% puv_pseudo_inv = pinv(puv)';
% 
% rxSignal = puv_pseudo_inv * EqualizedRxSymbolsMat;
% 
% % BER computation
% BER = comm.ErrorRate;
% 
% % Perform hard decision and measure errors
% rxBits = qamdemod(rxSignal, 2^bitsPerSubCarrier, 'OutputType', 'bit', ...
%     'UnitAveragePower', true);
% ber = BER(inpData(:), rxBits);
% 
% disp(['UFMC Reception, BER = ' num2str(ber(1)) ' at SNR = ' ...
%     num2str(snrdB) ' dB']);