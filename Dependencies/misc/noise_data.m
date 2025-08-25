function [DATA] = noise_data(DATA_unnoised,SNR)
%NOISE_DATA add gaussian distributed noise to the image, with defined SNR
%   SNR=mean(data)^2/std(noise)^2

if isempty(SNR)
    DATA=DATA_unnoised;
    fprintf('SNR is empty, no noise added \n')
else
    sig_power=mean(abs(DATA_unnoised).^2);
    noise_power=sig_power./SNR;
    
    real_noise=sqrt(noise_power/2)*randn(size(DATA_unnoised));
    imag_noise=sqrt(noise_power/2)*randn(size(DATA_unnoised)); %% complex
    %% noise the data
    DATA=DATA_unnoised + real_noise + 1i*imag_noise;
end
end

