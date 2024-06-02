clc
clear all
close all
set(0,'DefaultFigureWindowStyle','docked')

%% Set Up ability to save all the figures
save_fig = false;
if save_fig
    mkdir('./SimieleMattProjectPlots/')
end

%% Do Lowpass to Band Stop Filter Calculations from Report
% From Digital Filter Transformations Handout
fs = 1e6;
w_p = 0.2*pi;
wc_1 = (70000 * 2 * pi) / fs; wc_2 = (170000 * 2 * pi) / fs;
p = tan((wc_2 - wc_1) / 2) * tan(w_p/2);
lam = cos((wc_2 + wc_1) / 2) / cos((wc_2 - wc_1) / 2);

%% Define Analog Lowpass Filter

w_s1 = 2*atan((p*sin(0.18*pi)) / (cos(0.18*pi) - lam));
w_s2 = 2*atan((p*sin(0.24*pi)) / (cos(0.24*pi) - lam));
w_s = min(abs(w_s1), abs(w_s2));
% w_s = 0.400059*pi;
fs = 1e6;
Rp = 2; % Passband Ripple in dB, for type 2 this technically does not matter
Rs = 45; % Stopband Attenuation in dB
omega_p = tan(w_p/2); % Use an arbitrary Fs = 0.5 for the bilinear 
                      % transform
omega_s = tan(w_s/2); % Use an arbitrary Fs = 0.5 for the bilinear 
                      % transform

% Calculate Filter Order based on Chebyshev Equations
delt = 10^(2*Rp/20) - 1; 
A = 10^(Rs/20);

k = omega_p/omega_s;
k1 = delt/sqrt(A^2 - 1);

n = ceil(acosh(1/k1) / acosh(1/k));
% After some testing to get rid of the extra peak in the passband after 
% filter implementation a 5th order lowpass filter was used.

disp(['Low-Pass Filter Order Chosen: ', num2str(n)]); 
disp(['Band-Stop Filter Order: ', num2str(n*2)]); 
% Generate Analog Chebyshev Type II Filter 
[za,pa,ka] = cheby2(n,Rs,omega_s, 's');
plot_pz(save_fig, 's', za, pa, {'Analog Low Pass Filter', 'Pole/Zero Plot'})
[analog_num, analog_den] = zp2tf(za,pa,ka);
% Some versions of Matlab didn't like tf which is required to plot this
% plot_impulse(save_fig, 's', analog_num, analog_den, {'Analog Low Pass Filter', 'Filter Impulse Response'})

[h, w_out] = freqs(analog_num,analog_den,4096);
fig_fr = figure('Name', 'Analog Low Pass Frequency Response'); 
ax1 = subplot(2,1,1); semilogx(w_out, 20*log10(abs(h)), 'LineWidth', 2);
grid on;
xlabel('Frequency'); ylabel('Magnitude (dB)');
title({'Analog Low Pass Filter', 'Magnitude Response'});
set(gca, 'FontSize', 16)
ax2 = subplot(2,1,2); semilogx(w_out, angle(h), 'LineWidth', 2);
grid on;
xlabel('Frequency'); ylabel('Phase (radians)');
title({'Phase Response'});
set(gca, 'FontSize', 16)
linkaxes([ax1 ax2],'x')

if save_fig
    print(fig_fr,['./SimieleMattProjectPlots/', 'Analog Low Pass Frequency Response'],'-dpng')
end

%% Convert Analog Lowpass Filter to Digital Filter

[zdlp, pdlp, kdlp] = bilinear(za,pa,ka,0.5);
plot_pz(save_fig, 'z', zdlp, pdlp, {'Intermediate Digital Low Pass Filter', 'Pole/Zero Plot'})
[dlp_num, dlp_den] = zp2tf(zdlp,pdlp,kdlp);
plot_impulse(save_fig, 'z', dlp_num, dlp_den, {'Intermediate Digital Low Pass Filter', 'Impulse Response'})
plot_frequency_response(save_fig, dlp_num, dlp_den, {'Intermediate Digital Low Pass Filter', 'Frequency Response (using Fs)'}, fs)
plot_frequency_response(save_fig, dlp_num, dlp_den, {'Intermediate Digital Low Pass Filter', 'Frequency Response'})

%% Convert Digital Lowpass to Digital Band Stop

wp_hat = [0.14*pi, 0.34*pi];
ws_hat = [0.18*pi, 0.24*pi];
[dbs_num, dbs_den, transMatrix] = filter_z_transform('bs', dlp_num, dlp_den, lam, p);
% Transform to Pole Zero for For easy plotting and Possible 
% Better fixed point implementation
kdbs = dbs_num(1) / dbs_den(1);
k = kdbs;
zdbs = roots(dbs_num/dbs_num(1));
pdbs = roots(dbs_den/dbs_den(1));
plot_pz(save_fig, 'z', zdbs, pdbs, {'Floating Point Band-Stop Filter', 'Pole/Zero Plot'})
[dbsNum, dbsDen] = zp2tf(zdbs,pdbs,kdbs);
plot_impulse(save_fig, 'z', dbsNum, dbsDen, {'Floating Point Band-Stop Filter', 'Impulse Response'})
plot_frequency_response(save_fig, dbsNum, dbsDen, {'Floating Point Band-Stop Filter', 'Frequency Response (using Fs)'}, fs)
plot_frequency_response(save_fig, dbsNum, dbsDen, {'Floating Point Band-Stop Filter', 'Frequency Response'})

%% Setup Cascading filters
% Using a cascade of biquads here as it will help us be able to control the
% bounds of the input and output at specific locations in the filter
% structure
% 
% Get cascading filter coefficients where
%            A(z)       a(i,0) + a(i,1)z^-1 + a(i,2)z^-2
% H_i(z) = -------- = ------------------------------------
%            B(z)       b(i,0) + b(i,1)z^-1 + b(i,2)z^-2

phase_zeros = abs(angle(zdbs));
phase_poles = abs(angle(pdbs));
[~,indexPolesAngle] = sort(phase_poles,'descend');
[~,indexZerosAngle] = sort(phase_zeros,'descend');
zdbs = zdbs(indexZerosAngle);
pdbs = pdbs(indexPolesAngle);
N_filts = length(zdbs)/2;
num_cascade = zeros(N_filts, 3);
den_cascade = zeros(N_filts, 3);
for iFilt = 1:2:length(zdbs)
    num_cascade((iFilt+1)/2,:) = conv([1 -zdbs(iFilt)],[1 -zdbs(iFilt+1)]);
    den_cascade((iFilt+1)/2,:) = conv([1 -pdbs(iFilt)],[1 -pdbs(iFilt+1)]);
end

%% Convert to State Space Variables

A = zeros(2,2, size(num_cascade,1));
B = zeros(size(num_cascade,1),2);
C = zeros(size(num_cascade,1),2);
D = zeros(size(num_cascade,1),1);
for iFilt = 1:size(num_cascade,1)
    [A(:,:,iFilt), B(iFilt, :) , C(iFilt, :), D(iFilt)] = ...
        calc_state_space(num_cascade(iFilt, :), den_cascade(iFilt, :));
end

%% Plot the Fiter Design Floating Point Based on the State Matricies

[~, ~, filter_design_zeros, filter_design_poles, dbsFiltDesNum, dbsFiltDesDen] = ...
    get_filt_from_state_matrix(A, B, C, D, k);
plot_pz(save_fig, 'z', reshape(filter_design_zeros,[],1), reshape(filter_design_poles,[],1),...
    {'Floating Point Band-Stop Filter (State Matrix)', 'Pole/Zero Plot'})
plot_impulse(save_fig, 'z', dbsFiltDesNum, dbsFiltDesDen, {'Floating Point Band-Stop Filter (State Matrix)', 'Impulse Response'})
plot_frequency_response(save_fig, dbsFiltDesNum, dbsFiltDesDen, {'Floating Point Band-Stop Filter (State Matrix)', 'Frequency Response (Using Fs)'}, fs)
plot_frequency_response(save_fig, dbsFiltDesNum, dbsFiltDesDen, {'Floating Point Band-Stop Filter (State Matrix)', 'Frequency Response'})

%% Calculate SNR and Required number of bits

[noise_out, gain_per_ss] = compute_snr(A,B,C,D,k);
% Sort from higher noise variance in State Matricies
% Should yield a lower noise power level since quantization noise is not
% compounded over the cascade
[~,indexGainDescend] = sort(sum(gain_per_ss,2),'descend');
A(:,:,:) = A(:,:,indexGainDescend);
B(:,:) = B(indexGainDescend,:);
C(:,:) = C(indexGainDescend,:);
[noise_out_des, gain_per_ss_des] = compute_snr(A,B,C,D,k);
% D is always one in this case so leave it alone here
L = ceil( (log2(noise_out * 10^(7.2)) + 3) / 2 );

T = array2table(horzcat(reshape(A,[],size(A,3),1)', B,C,D),...
    'VariableNames',{'a11','a12','a21', 'a22', 'b1', 'b2', 'c1', 'c2', 'd'});

disp('----------------------------------------------------------------------------');
disp('Filter Coefficients Floating Point')
disp(T)
disp('----------------------------------------------------------------------------');

disp(['Number of Bits Used in Fixed Point Representation: ', num2str(L+1)]);

%% Convert the state matricies to fixed point represetation
% Note float used here is the decimal representation of the fixed point
% number. This is limited to L+1 bits where the +1 bit represents the signs
[A_fixed, A_float] = convert_to_fixed(A,L);
[B_fixed, B_float] = convert_to_fixed(B,L);
[C_fixed, C_float] = convert_to_fixed(C,L);
% Since D is 1 for all of the cascades this does not need to be converted
% to fixed there will technically be no "multiplier" d
[k_fixed, k_float] = convert_to_fixed(k,L);

T_fixed = array2table(horzcat(reshape(A_float,[],size(A,3),1)', B_float,C_float,D),...
    'VariableNames',{'a11','a12','a21', 'a22', 'b1', 'b2', 'c1', 'c2', 'd'});

disp('----------------------------------------------------------------------------');
disp('Filter Coefficients Fixed Point')
disp(T_fixed)
disp('----------------------------------------------------------------------------');

%% Get the Quantized Filter

[quantized_num, quantized_den, quantized_zeros_filt, quantized_poles_filt, total_num_quant, total_den_quant] = ...
    get_filt_from_state_matrix(A_float, B_float, C_float, D, k_float);
plot_pz(save_fig, 'z', reshape(quantized_zeros_filt,[],1), reshape(quantized_poles_filt,[],1),...
    {'Fixed Point Band-Stop Filter', 'Pole/Zero Plot'})
plot_impulse(save_fig, 'z', total_num_quant, total_den_quant, {'Fixed Point Band-Stop Filter', 'Impulse Response'})
plot_frequency_response(save_fig, total_num_quant, total_den_quant, {'Fixed Point Band-Stop Filter', 'Frequency Response (using Fs)'}, fs)
plot_frequency_response(save_fig, total_num_quant, total_den_quant, {'Fixed Point Band-Stop Filter', 'Frequency Response'})

%% Time Domain SNR Simulation
% Generate a sinusoid of -6dB full scale within pass band of the filter
%L = L-1;
t = 0:0.01:15;
x = 0.5*cos(0.15*pi*t);
[y, y_1, y_2, y_3]= ...
    time_domain_sim(A,B,C,D,k,x);
[y_quant, y_1_quant, y_2_quant, y_3_quant] = ...
    time_domain_sim(A_float, B_float, C_float, D, k_float, x, L);
simulated_snr = 10*log10( sum(y(end,3:end).^2) / sum(abs(y(end,3:end) - y_quant(end,3:end)).^2) );
disp(['Simulated SNR (dB): ', num2str(simulated_snr)]);

figure('Name', 'Time Domain Simulation Results')
yyaxis left
plot(t, y(length(D)+1, 3:end)); hold on;
plot(t, y_quant(length(D)+1, 3:end), 'r'); hold off;
legend('Floating Point', 'Fixed Point')
xlabel('Samples');
ylabel('Output Amplitude');
ylim([-1 1])
yyaxis right
plot(t, abs(y_quant(length(D)+1, 3:end) - y(length(D)+1, 3:end)));
ylabel('Sample Error in Amplitude')
ylim([0 0.001])
legend('Floating Point', 'Fixed Point', 'Difference')
title('Time Domain Analysis')
grid on;
set(gca, 'FontSize', 16)

% Compute the Used Adders and Mutipliers Etc:
[total_basic_logic] = calc_hardware_usage(1, 0, length(D), 3, 8, 2, L);
disp(['Basic Logic Gates Used: ', num2str(total_basic_logic)]);

function [y, y_1, y_2, y_3] = time_domain_sim(A,B,C,D,k,x,L)
    % Store output of each filter in a cascade with y(i,:) being the inital 
    % input to the filter structure
    y = zeros(length(D)+1, length(x)+2);
    y(1,3:end) = x;
    % Store output of each State Matrix and the "d" path in a cascade 
    y_1 = zeros(length(D)+1, length(x)+2); % "d" path
    y_2 = zeros(length(D)+1, length(x)+2); % State Matrix 1 (eq u_1 in report)
    y_3 = zeros(length(D)+1, length(x)+2); % State Matrix 2 (eq u_2 in report)

    % If bit length is not defined then use floating point numbers
    if nargin == 6
        % Scale x by input gain
        y(1,:) = y(1,:)*k;
        % Start calulating output when we have enough inputs 
        for sampIdx = 3:(length(x)+2)
            for filtIdx = 1:length(D)
                y_1(filtIdx+1,sampIdx) = y(filtIdx,sampIdx)*D(filtIdx);
                y_2(filtIdx+1,sampIdx) = ((B(filtIdx,1)*y(filtIdx,sampIdx-1)) + ...
                    ((A(1,2,filtIdx)*B(filtIdx,2) - A(2,2,filtIdx)*B(filtIdx,1)) * y(filtIdx,sampIdx-2))) - ...
                    (-1* trace(A(:,:,filtIdx))*y_2(filtIdx+1,sampIdx-1) + det(A(:,:,filtIdx))*y_2(filtIdx+1,sampIdx-2));
                y_3(filtIdx+1,sampIdx) = ((B(filtIdx,2)*y(filtIdx,sampIdx-1)) + ...
                    ((A(2,1,filtIdx)*B(filtIdx,1) - A(1,1,filtIdx)*B(filtIdx,2)) * y(filtIdx,sampIdx-2))) - ...
                    (-1* trace(A(:,:,filtIdx))*y_3(filtIdx+1,sampIdx-1) + det(A(:,:,filtIdx))*y_3(filtIdx+1,sampIdx-2));
                y(filtIdx+1, sampIdx) = y_1(filtIdx+1,sampIdx) + ...
                    C(filtIdx,1) * y_2(filtIdx+1,sampIdx) + ...
                    C(filtIdx,2) * y_3(filtIdx+1,sampIdx);
            end
        end
    else 
        % Scale x by input gain
        y(1,:) = Q_comp(y(1,:)*k,L);
        % Start calulating output when we have enough inputs 
        for sampIdx = 3:(length(x)+2)
            for filtIdx = 1:length(D)
                y_1(filtIdx+1,sampIdx) = Q_comp(y(filtIdx,sampIdx)*D(filtIdx),L);
                y_2(filtIdx+1,sampIdx) = Q_comp((B(filtIdx,1)*y(filtIdx,sampIdx-1)),L)...
                    + Q_comp(Q_comp((A(1,2,filtIdx)*B(filtIdx,2)),L) * y(filtIdx,sampIdx-2),L)...
                    - Q_comp(Q_comp((A(2,2,filtIdx)*B(filtIdx,1)),L) * y(filtIdx,sampIdx-2),L)...
                    + Q_comp(A(1,1,filtIdx) * y_2(filtIdx+1,sampIdx-1),L)...
                    + Q_comp(A(2,2,filtIdx) * y_2(filtIdx+1,sampIdx-1),L)...
                    - Q_comp(Q_comp(A(2,2,filtIdx)*A(1,1,filtIdx),L) * y_2(filtIdx+1,sampIdx-2),L)...
                    + Q_comp(Q_comp(A(1,2,filtIdx)*A(2,1,filtIdx),L) * y_2(filtIdx+1,sampIdx-2),L);
                y_2(filtIdx+1, sampIdx) = Q_comp(y_2(filtIdx+1, sampIdx),L);
                y_3(filtIdx+1,sampIdx) = Q_comp((B(filtIdx,2)*y(filtIdx,sampIdx-1)),L) + ...
                    Q_comp(Q_comp((A(2,1,filtIdx)*B(filtIdx,1)),L) * y(filtIdx,sampIdx-2),L) - ...
                    Q_comp(Q_comp(A(1,1,filtIdx)*B(filtIdx,2),L) * y(filtIdx,sampIdx-2),L) + ...
                    Q_comp(A(1,1,filtIdx) * y_3(filtIdx+1,sampIdx-1),L) + ...
                    Q_comp(A(2,2,filtIdx) * y_3(filtIdx+1,sampIdx-1),L) - ...
                    Q_comp(Q_comp(A(2,2,filtIdx)*A(1,1,filtIdx),L) * y_3(filtIdx+1,sampIdx-2),L) + ...
                    Q_comp(Q_comp(A(1,2,filtIdx)*A(2,1,filtIdx),L) * y_3(filtIdx+1,sampIdx-2),L);
                y_3(filtIdx+1, sampIdx) = Q_comp(y_3(filtIdx+1, sampIdx),L);
                y(filtIdx+1, sampIdx) = Q_comp(y_1(filtIdx+1,sampIdx),L) + ...
                    Q_comp(C(filtIdx,1) * y_2(filtIdx+1,sampIdx),L) + ...
                    Q_comp(C(filtIdx,2) * y_3(filtIdx+1,sampIdx),L);
                y(filtIdx+1, sampIdx) = Q_comp(y(filtIdx+1, sampIdx),L);
            end
        end
    end
end

%% Hardware Calculation

function [total_basic_logic] = calc_hardware_usage(input_mult, output_mult, ...
    num_cascades, num_add_cas, num_mult_cas, num_reg_cas, num_bits_mag)

    total_bits = num_bits_mag + 1;

    total_basic_logic = (input_mult + output_mult) * 2 * total_bits^2;

    total_basic_logic = total_basic_logic ...
        + (num_cascades * num_add_cas * 4 * total_bits) ...
        + (num_cascades * num_mult_cas * 2 * total_bits^2) ...
        + (num_cascades * num_reg_cas * 5 * total_bits);

end

%% Additional Functions for Project Script

% - Start -- Functions Used for Fixed Point Implementation ------------- %
function [fixed_point_arr, floating_point_arr] = convert_to_fixed(array2quantize, mag_length)
% set scaling factor to scale the number to fixed length 
scale_factor = 2^mag_length;
% Get Max value
max_magnitude = 2^mag_length - 1;

fixed_point_arr = zeros(size(array2quantize));
% Convert to a fixed point represetation
for i = 1:size(array2quantize,1)
    for j = 1:size(array2quantize,2)
        for k = 1:size(array2quantize,3)
            scaled_num = abs(array2quantize(i,j,k)) * scale_factor;
            rounded_num = round(scaled_num);
            clamped_magnitude = min(rounded_num, max_magnitude);
            if rounded_num > max_magnitude 
                disp('Rounding Error incurred')
            end
            % Determine the sign and apply sign-magnitude representation
            if array2quantize(i,j,k) < 0
                sign_bit = 2^mag_length; % Set the sign bit (equivalent to adding 1 in the MSB)
                fixed_point_arr(i,j,k) = sign_bit + clamped_magnitude;
            else
                fixed_point_arr(i,j,k) = clamped_magnitude;
            end
        end
    end
end
floating_point_arr = zeros(size(array2quantize));
for i = 1:size(fixed_point_arr, 1)
    for j = 1:size(fixed_point_arr, 2)
        for k = 1:size(fixed_point_arr, 3)
            fixed_value = fixed_point_arr(i,j,k);
            if bitget(fixed_value, mag_length + 1) % Check the sign bit
                % Negative number
                magnitude = bitset(fixed_value, mag_length + 1, 0); % Clear the sign bit
                floating_point_arr(i,j,k) = -magnitude / scale_factor;
            else
                % Positive number
                floating_point_arr(i,j,k) = fixed_value / scale_factor;
            end
        end
    end
end

end

function quant_out = Q_comp(x,L)
    [~, quant_out] = convert_to_fixed(x,L);
end
% -- End --- Functions Used for Fixed Point Implementation ------------- %

% State Space Transformations SNR and Transforming back into a filter
function [a, b, c, d] = calc_state_space(num, den)
a = zeros(2,2);
b = zeros(1,2);
c = zeros(1,2);

% Use Equations from given papers
d = num(1);
a(1:size(a,1)+1:end) = -den(2)/2;
a(2,1) = sqrt( den(3) - ((den(2)^2) / 4) );
a(1,2) = -a(2,1);

alpha = num(2) - (d * den(2));
beta = num(3) - (d * den(3));

gamma = ((beta + (alpha*a(2,2)))^2) / (4*a(1,2)^2);
b(2) = sqrt( ( sqrt(alpha^2 + (4 * gamma ) ) - alpha ) / 2); c(2) = -1*b(2);
c(1) = -1 * sqrt(alpha + c(2)^2); b(1) = c(1);
end

function [noise_out, gain_per_state_space] = compute_snr(A,B,C,D,k)
M = length(D);
noise_out = k^2;
gain_per_state_space = zeros(length(D),2);
for i = 1:M
    state_space_1 = [0, B(i,1), (A(1,2,i)*B(i,2) - A(2,2,i)*B(i,1));...
        1, -1*trace(A(:,:,i)), det(A(:,:,i))];
    
    state_space_2 = [0, B(i,2), (A(2,1,i)*B(i,1) - A(1,1,i)*B(i,2));...
        1, -1*trace(A(:,:,2)), det(A(:,:,i))];
    
    H_num = [1];
    H_den = [1];
    for j = (i+1):M
        A_0 = D(j);
        A_1 = B(j,1)*C(j,1) + B(j,2)*C(j,2) - D(j)*trace(A(:,:,j));
        A_2 = D(j)*det(A(:,:,j)) + C(j,2)*B(j,1)*A(2,1,j) + C(j,1)*B(j,2)*A(1,2,j)...
            - C(j,1)*B(j,1)*A(2,2,j)-C(j,2)*B(j,2)*A(1,1,j);
        B_1 = -1*trace(A(:,:,j));
        B_2 = det(A(:,:,j));
        den = [1, B_1, B_2];
        num = [A_0, A_1, A_2];
        H_num = conv(H_num, num);
        H_den = conv(H_den, den);
    end
    gain_per_state_space(i,1) = filternorm(state_space_1(1,:), state_space_1(2,:),2);
    gain_per_state_space(i,2) = filternorm(state_space_2(1,:), state_space_2(2,:),2);
    % This equation comes from More economical state-space digital-filter structures 
    % which are free of constant-input limit cycles
    % 2 sources of Q Error are input to the State Matrices and 3 more add
    % into the full transfer function (could make this 2 in my case since d
    % is never a "multiplier"
    noise_out = noise_out + (2*filternorm(state_space_1(1,:), state_space_1(2,:),2)^2 + 2*filternorm(state_space_2(1,:), state_space_2(2,:),2)^2+ 3) ...
        * filternorm(H_num, H_den,2)^2; 
end
% Multiply by q^2/12 but just 1/12 here as we know q^2 is in the sinusoid
% used for the SNR calculation
noise_out = noise_out/12;
end

function [num, den, zeros_filt, poles_filt, total_num, total_den] = get_filt_from_state_matrix(A,B,C,D,k)
    M = length(D);
    num = zeros(length(D),3);
    den = zeros(length(D),3);
    zeros_filt = zeros(length(D),2);
    poles_filt = zeros(length(D),2);
    total_num = [1];
    total_den = [1];
    for j = 1:M
        A_0 = D(j);
        A_1 = B(j,1)*C(j,1) + B(j,2)*C(j,2) - D(j)*trace(A(:,:,j));
        A_2 = D(j)*det(A(:,:,j)) + (C(j,2)*B(j,1)*A(2,1,j)) + (C(j,1)*B(j,2)*A(1,2,j)) ...
            - (C(j,1)*B(j,1)*A(2,2,j)) - (C(j,2)*B(j,2)*A(1,1,j));
        B_1 = -1*trace(A(:,:,j));
        B_2 = det(A(:,:,j));
        den(j,:) = [1, B_1, B_2];
        num(j,:) = [A_0, A_1, A_2];
        zeros_filt(j,:) = roots(num(j,:))';
        poles_filt(j,:) = roots(den(j,:))';
        total_num = conv(total_num, num(j,:));
        total_den = conv(total_den, den(j,:));
    end
    total_num = k*total_num;
end

% ----- Wrappers for Plotting Functions used in the script ------------- %
function plot_pz(save_fig, domain, zeros, poles, plot_title)

    fig_pz = figure('Name', string(plot_title{1}) + ' ' + string(plot_title{2}));
    scatter(real(zeros), imag(zeros), 100,'filled', 'g', 'o'); hold on;
    scatter(real(poles), imag(poles), 100, 'r', 'x'); hold off;
    xlabel('real'); ylabel('imaginary');
    minx = floor(min(min(real(zeros)), min(real(poles))));
    miny = floor(min(min(imag(zeros)), min(imag(poles))));
    maxx = ceil(max(max(real(zeros)), max(real(poles))));
    maxy = ceil(max(max(imag(zeros)), max(imag(poles))));
    grid on;
    title(plot_title)

    if domain == 's'
        hold on;
        plot([0,0], [miny, maxy], '--k');
        hold off;
        xlim([minx, maxx]);
        ylim([miny, maxy]);
    end
    if domain == 'z'
        hold on;
        th = 0:pi/50:2*pi;
        xunit = cos(th);
        yunit = sin(th);
        plot(xunit, yunit, '--k');
        hold off;
        xlim([-1.5, 1.5]);
        ylim([-1.5, 1.5]);
    end
    legend('zeros', 'poles')
    set(gca, 'FontSize', 16)
    if save_fig
        plot_save = strrep(string(plot_title{1}) + ' ' + string(plot_title{2}),'/','-');
        plot_save = "./SimieleMattProjectPlots/" + plot_save;
        print(fig_pz,plot_save,'-dpng')
    end
end

function plot_impulse(save_fig, domain, num, den, plot_title)

    fig_impulse = figure('Name', string(plot_title{1}) + ' ' + string(plot_title{2}));

    if domain == 's'
        sys = tf(num, den);
        [y, tout] = impulse(sys);
        stairs(tout, y, 'LineWidth', 2)
        xlabel('time (seconds)');
        ylabel('Amplitude');
    end
    if domain == 'z'
        [y, tout] = impz(num, den);
        stem(tout, y, 'LineWidth', 2)
        xlabel('Samples');
        ylabel('Amplitude');
    end
    title(plot_title)
    set(gca, 'FontSize', 16)
    if max(tout) > 200
        xlim([0,100])
    end
    
    if save_fig
        plot_save = strrep(string(plot_title{1}) + ' ' + string(plot_title{2}),'/','-');
        plot_save = "./SimieleMattProjectPlots/" + plot_save;
        print(fig_impulse,plot_save,'-dpng')
    end
end

function plot_frequency_response(save_fig, num, den, plot_title, fs)
    if nargin == 4
        fig_fr = figure('Name', string(plot_title{1}) + ' ' + string(plot_title{2}));
        [h,w] = freqz(num, den);
        ax1 = subplot(2,1,1); plot(w/pi, 20*log10(abs(h)), 'LineWidth', 2);
        grid on;
        xlabel('Normalized Frequency (\times \pi rad/sec)'); ylabel('Magnitude (dB)');
        title({string(plot_title{1}), 'Magnitude Response'});
        set(gca, 'FontSize', 16)
        ax2 = subplot(2,1,2); plot(w/pi, angle(h), 'LineWidth', 2);
        grid on;
        xlabel('Normalized Frequency (\times \pi rad/sec)'); ylabel('Phase (radians)');
        title({'Phase Response'});
        set(gca, 'FontSize', 16)
        linkaxes([ax1 ax2],'x')
    elseif nargin == 5
        fig_fr = figure('Name', string(plot_title{1}) + ' ' + string(plot_title{2}));
        [h,w] = freqz(num, den, 16384, fs);
        ax1 = subplot(2,1,1); plot(w/1e3, 20*log10(abs(h)), 'LineWidth', 2);
        grid on;
        xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
        title({string(plot_title{1}), 'Magnitude Response'});
        set(gca, 'FontSize', 16)
        ax2 = subplot(2,1,2); plot(w/1e3, angle(h), 'LineWidth', 2);
        grid on;
        xlabel('Frequency (kHz)'); ylabel('Phase (radians)');
        title({'Phase Response'});
        set(gca, 'FontSize', 16)
        linkaxes([ax1 ax2],'x')
    end
    if save_fig
        plot_save = strrep(string(plot_title{1}) + ' ' + string(plot_title{2}),'/','-');
        plot_save = "./SimieleMattProjectPlots/" + plot_save;
        print(fig_fr,plot_save,'-dpng')
    end
end
%% Function to do Digital Transformations

function [num_out, den_out, transMatrix] = filter_z_transform(filter_type, num_in, den_in, c1, c2)
    % Computes the numerator, denominator, and transformation matrix for a given 
    % digital filter transformation 
    
    % Check if num_in and den_in have the same order
    if length(num_in) ~= length(den_in)
        error('Transfer function does not have the same order in den and numerator!');
    else 
        order = length(num_in)-1;
    end

    if strcmp(filter_type, 'lp')
        [num_trans, den_trans] = lp2lp_z_transform(c1);
    elseif strcmp(filter_type, 'hp')
        [num_trans, den_trans] = lp2hp_z_transform(c1);
    elseif strcmp(filter_type, 'bp')
        [num_trans, den_trans] = lp2bp_z_transform(c1, c2);
    elseif strcmp(filter_type, 'bs')
        [num_trans, den_trans] = lp2bs_z_transform(c1, c2);
    else
        error('Desired filter type is not supported. Please use: lp, hp, bp, or bs!');
    end

    transMatrix = compute_transfer_matrix(num_trans, den_trans, order);
    
    [num_out, den_out] = compute_output_transfer_func(transMatrix, num_in, den_in);
end

function [num_trans, den_trans] = lp2lp_z_transform(c1)
    num_trans = [-1*c1, 1];
    den_trans = [1, -1*c1];
end

function [num_trans, den_trans] = lp2hp_z_transform(c1)
    num_trans = [-1*c1, -1];
    den_trans = [1, c1];
end

function [num_trans, den_trans] = lp2bp_z_transform(c1, c2)
    num_trans = -1.*[(c2-1), -2*c1*c2, (c2+1)];
    den_trans = [(c2+1), -2*c1*c2, (c2-1)];
end

function [num_trans, den_trans] = lp2bs_z_transform(c1, c2)
    num_trans = [(1-c2), -2*c1, (1+c2)];
    den_trans = [(1+c2), -2*c1, (1-c2)];
end

function transMatrix = compute_transfer_matrix(num, den, order)
    transMatrix = zeros(order+1, (order * (length(num) - 1)) + 1);
    for i = 1:order+1
        num_conv = 1;
        den_conv = 1;
        for j = 1:i-1
            num_conv = conv(num_conv, num);
        end
        for j = 1:order - i + 1
            den_conv = conv(den_conv, den);
        end
        transMatrix(i,:) = conv(num_conv, den_conv);
    end
end

function [num_out, den_out] = compute_output_transfer_func(transMatrix, num, den)
    num_out = zeros(size(transMatrix(1,:)));
    den_out = zeros(size(transMatrix(1,:)));
    for i = 1:length(num)
       den_out = den_out + den(i) * transMatrix(i,:);
       num_out = num_out + num(i) * transMatrix(i,:);  
    end
    num_out = num_out / den_out(1);
    den_out = den_out / den_out(1);
end
