function result = radial_filter_fft(input, fft_kernel,  km, kn)
    % Fast convolution using FFT
    [m, n] = size(input);
 

    % Zero-pad kernel and input to same size
    pad_m = m + km - 1;
    pad_n = n + kn - 1;

    % FFT-based convolution
    fft_input = fft2(input, pad_m, pad_n);
    result_full = real(ifft2(fft_input .* fft_kernel));

    % Crop to original size
    start_row = floor(km/2);
    start_col = floor(kn/2);
    result = result_full(start_row+(1:m), start_col+(1:n));
end