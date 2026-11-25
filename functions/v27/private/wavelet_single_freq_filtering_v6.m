function [filtered_data, filtered_complex] = wavelet_single_freq_filtering_v6(data, srate, freq, numcycles)
% Custom Morlet convolution; preserve the legacy numerical calculation.
% Restore the caller warning settings after temporary legacy suppression.
    previousWarnings = warning;
    restoreWarnings = onCleanup(@() warning(previousWarnings)); %#ok<NASGU>
    warning('off', 'all');

    center_freq = freq;

    fwhm = numcycles * sqrt(2 * log(2)) / (pi * center_freq);
    time = -3 * fwhm : 1 / srate : 3 * fwhm;

    wavelet = exp(2 * 1i * pi * center_freq .* time) .* exp(-4 * log(2) * time.^2 / fwhm^2);
    wavelet = wavelet / sum(abs(wavelet));

    n_wavelet = length(time);
    n_data = size(data, 1);
    n_ch = size(data, 2);
    n_convolution = 2^nextpow2(n_wavelet + n_data - 1);
    half_of_wavelet_size = (n_wavelet - 1) / 2;

    fft_wavelet = fft(wavelet, n_convolution);
    fft_data = fft(data', n_convolution);
    convolution_result_fft = ifft(fft_wavelet .* fft_data, n_convolution, 2);
    convolution_result_fft = convolution_result_fft(:, half_of_wavelet_size+1 : half_of_wavelet_size+n_data)';

    filtered_complex = convolution_result_fft;
    filtered_data = real(convolution_result_fft);
end
