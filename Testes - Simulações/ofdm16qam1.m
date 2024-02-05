% Eb = linspace(5, 15, 11);
% BER = [];
% for i = 1:size(Eb, 2)
%     EbNo = Eb(i);
%     sim('ofdm.slx');
%     BER(i) = simout(1);
% end
% 
% semilogy(Eb, BER);
% % plot(Eb, BER);
% % title('OFDM 16-QAM BER vs Eb/No')
% % xlabel('Eb/No')
% % ylabel('BER')

snr = [0 3 6 9 12];
ber_16qam_eb_no = [0.3854 0.3854 0.3646 0.2135 0.1667];

semilogy(snr, ber_16qam_eb_no)