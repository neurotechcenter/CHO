function [ filtered_data filtered_complex ] = wavelet_single_freq_filtering_v4(data, srate, freq, numcycles)
% Custom Morlet convolution; preserve the legacy numerical calculation.
% Restore the caller warning settings after temporary legacy suppression.


center_freq = freq;
w_sigma = numcycles/(2*pi*center_freq);
time        = -w_sigma*3:1/srate:3*w_sigma;


fwhm = numcycles * (2*log(2)).^0.5/(pi*center_freq);
wavelet = exp(2*1i*pi*center_freq.*time).*exp(-4*log(2)*time.^2/fwhm.^2);

wavelet = wavelet / sum(abs(wavelet));
previousWarnings = warning;
restoreWarnings = onCleanup(@() warning(previousWarnings)); %#ok<NASGU>
warning('off', 'all');
half_of_wavelet_size = ((length(time)-1)/2);

n_wavelet     = length(time);
n_data        = size(data,1);
n_ch = size(data,2);
n_convolution = n_wavelet+n_data-1;

fft_wavelet = fft(wavelet,n_convolution);

filtered_data  = zeros(n_data,n_ch);
filtered_complex = zeros(n_data,n_ch);

for chani=1:n_ch
    fft_data = fft(squeeze(data(:,chani))',n_convolution);
    convolution_result_fft = ifft(fft_wavelet.*fft_data,n_convolution) * sqrt(4/(2*pi*center_freq));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);    
    filtered_complex(:,chani)  = convolution_result_fft;
    filtered_data(:,chani)  = real(convolution_result_fft);

end


