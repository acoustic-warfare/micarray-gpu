function audio_out = AW_listening_parallel(c,f_sampling,r_prime,theta_listen,phi_listen,audio_signals,filter_coefficients,center_frequencies,weight_m)
%
%   Spatial filtering in r_hat(theta,phi)
%
%
%   This function takes in the audio_signals generated and the
%   corresponding array matrix-geometries. This function further needs some
%   information about the audio_signals such as the sampling frequency,
%   soundspeed in the medium (c).
%
%   A direction is needed (theta,phi) to output a spatially filtered
%   signal.
%
%   The function splits the input audio signals recorded by the
%   array_matrices and then perform narrowband beamforming on the signals.
%   Different audio signals corresponds to different sub_arrays.
%

% Initialize our audio out vector
samples = length(audio_signals(:,1));
audio_out = zeros(samples,1);

% Convert input theta and phi in degrees to radians
% theta = pi/180*theta_listen;
% phi = pi/180*phi_listen;

uni_distance = abs(r_prime(2,1) - r_prime(2,2));

theta = theta_listen;
phi = phi_listen;

x_factor = sin(theta)*cos(phi);
y_factor = sin(theta)*sin(phi);

% Filter the signals generated on the arrays
x = r_prime(1,:);
y = r_prime(2,:);
elements = length(r_prime(1,:));

for freq_ind = 1:length(filter_coefficients(:,1))
    % Filter coefficients for each band
    b = filter_coefficients(freq_ind,:);

    % Center frequency
    frequency = center_frequencies(freq_ind);

    % The narrowband frequency
    k = 2*pi*frequency/c;

    % The normalized frequency
    ny = frequency/f_sampling;

    w_index = weight_index_parallel(uni_distance,frequency,c);
    weight = weight_m(w_index,:);

    audio_temp = zeros(samples,1); 
    mic_data = zeros(samples,1);

    for (mic_ind = 1:elements)
        if weight(mic_ind) == 1
            % Apply filter on every signal recorded from the elements
            audio_temp(:,1) = filter(b,1,audio_signals(:,mic_ind));

            % Phase shift value that is dependent on the frequency and
            % the location of the element (x,y)
            phase_shift_value = -k*(x(mic_ind)*x_factor + ...
                y(mic_ind)*y_factor);
            
            % Add the signals together
            mic_data(:,1) = mic_data(:,1) + ...
                phase_shift(audio_temp,ny,phase_shift_value);

            %mic_data(:,mic_ind) = phase_shift(filter(b,1,audio_signal(:,mic_ind)),ny,phase_shift_value);
        end
    end
   % audio_out = audio_out + sum(mic_data,2);
   audio_out = audio_out + mic_data;
end
    
%audio_out = audio_out/max(audio_out);
end
