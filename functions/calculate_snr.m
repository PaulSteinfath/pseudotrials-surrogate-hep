function snr_stats = calculate_snr(original_data, p300_amplitude, hep_amplitude)
    % CALCULATE_SNR Computes the Signal-to-Noise Ratio (SNR) for P300 and HEP signals.
    %
    %   SNR_STATS = CALCULATE_SNR(ORIGINAL_DATA, P300_AMPLITUDE, HEP_AMPLITUDE)
    %   calculates the SNR for P300 and HEP signals based on the provided
    %   original data and their respective amplitudes. The function returns a
    %   structure containing the SNR values in decibels (dB) for both P300 and
    %   HEP signals.
    %
    %   Inputs:
    %       ORIGINAL_DATA - The original data array containing the signal and noise.
    %       P300_AMPLITUDE - The amplitude values of the P300 signal.
    %       HEP_AMPLITUDE - The amplitude values of the HEP signal.
    %
    %   Outputs:
    %       SNR_STATS - A structure with fields:
    %           'snr_p300' - SNR of the P300 signal in dB.
    %           'snr_hep'  - SNR of the HEP signal in dB.

    % Calculate signal and noise power
    signal_power_p300 = mean(p300_amplitude.^2);
    signal_power_HEP = mean(hep_amplitude.^2);
    noise_power = mean(reshape(original_data, 1, []).^2);
    
    % Calculate SNR in dB
    snr_stats = struct(...
        'snr_p300', 10 * log10(signal_power_p300 / noise_power), ...
        'snr_hep', 10 * log10(signal_power_HEP / noise_power) ...
    );
end