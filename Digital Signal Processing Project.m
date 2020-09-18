clear;
clc;

%import of the audio track
[y, fs] = audioread('Input Audio File.wav');
N = length(y);
d = N/fs;
fprintf('Sampling frequency = %d Hz \n', fs);
fprintf('Track duration = %5.2f s \n\n', d);

%we plot the spectrum before the computation
Y = fft(y(1:N));
Y_norm = Y/N;
Y_dB = 20*log10(abs(Y));
f = linspace(0, fs, N);
t = linspace(1/fs, d, N);

figure(1)
subplot(2,1,1);              
plot(t, y, 'Color', [0.6350, 0.0780, 0.1840]);                                     
title('Original signal: time domain'); %time domain
xlabel('nT (s)'); ylabel('y(nT)');
axis([0.1 d min(y) max(y)])
grid;
subplot(2,1,2);                                        
plot(f, Y_dB, 'Color', [0.6350, 0.0780, 0.1840]); 
title('Original signal: DFT'); %frequency domain
xlabel('f (Hz)'); ylabel('|Y(f)| (dB)');
maxy = max(Y_dB); miny = maxy - 200;
axis([0 fs miny maxy]);
grid;

%Now we find the frequencies f_1 and f_2 and their respective amplitudes 
%A1 and A2 searching peaks
peaks = find(abs(Y) >= (max(abs(Y))/2));
A1 = abs(Y_norm(peaks(1))) + abs(Y_norm(peaks(4)));
A2 = abs(Y_norm(peaks(2))) + abs(Y_norm(peaks(3)));

freq = peaks * fs / N;

fprintf('Estimate of the frequencies and amplitudes of the carriers: \n\n');
fprintf('f_1 = %d Hz \n', freq(1));
fprintf('f_2 = %d Hz \n\n', freq(2));
fprintf('A_1 = %5.4f\n', A1);
fprintf('A_2 = %5.4f\n', A2);

%Parameters of the passband filter: as the assignment suggest, we'll 
%implement a general IIR filter with N = 2, M = 0, finding manually the
%parameters b0, a1 and a2 (as we've done during the course)
%We first define the general parameters of filters, taking a very narrow
%band
f_3dB = 1;                 
theta_3dB = 2*pi*f_3dB/fs; 

%As we know, for a passband: theta = 2*delta => delta = theta/2,
%2*delta = 2*(1 - r) => r = 1 - delta, and we can approximate
%b0 as delta = (1 - r) (so we normalize it to 1). In this case we're using
%the improvement filter design we've studied during the lessons.
delta = theta_3dB / 2;
r = 1 - delta;

%Now we need to find the parameters of our filter with N = M = 2: 
%from the theory, setting p1 and p2 as complex conjugate poles, we have 
%that a1 = 2*r*(cos(theta_zero)) and a2 = - r*r; then we extract the two
%carriers. At the very end of the process, we compute the DFT and the 
%its normalization 
a1 = zeros(2);
a2 = zeros(2);
b0 = zeros(2);

for i = 1:2
    a1(i) = 2*r*cos(2*pi*freq(i)/fs);
    a2(i) = -r*r;
    b0(i) = delta*sin(2*pi*freq(i)/fs);
end

%carrier one
num_1 = b0(i);
den_1 = [1 -a1(1) -a2(1)];
[H_1, w] = freqz(num_1, den_1, 'whole', 2048, fs);
carry1 = filter(num_1, den_1, y);
C1 = fft(carry1(1:N)); 
C1_norm = C1 / N;
C1_dB = 20*log10(abs(C1_norm));

%carrier two
num_2 = b0(2);
den_2 = [1 -a1(2) -a2(2)];
[H_2, w] = freqz(num_2, den_2, 'whole', 2048, fs);
carry2 = filter(num_2, den_2, y);
C2 = fft(carry2(1:N)); 
C2_norm = C2 / N;
C2_dB = 20*log10(abs(C2_norm));

%Plot of the spectrum of the two carriers
figure(2)                                   
plot(f, C1_dB, f, C2_dB);
title('Magnitude (dB) of the spectrum of carriers');
xlabel(' f (Hz)');
legend({'Carrier one','Carrier two'},'Location','northeast')
maxy = max(Y_dB); miny = maxy - 350;
axis([0 fs miny maxy]);
grid;

%Plot of the filters used to extract carriers
figure(3)                       
f_2048 = linspace(0, fs, 2048);          
subplot(2,2,1);
plot(f_2048, 20*log10(abs(H_1)))
xlabel(' f (Hz)'); ylabel('|H_1(f)|  (dB)');
maxy = max(20*log10(abs(H_1))) + 10; miny = maxy - 200;
axis([0 fs miny maxy]);
grid;
subplot(2,2,2);
plot(f_2048, 20*log10(abs(H_2)), 'Color', [0.8500, 0.3250, 0.0980])
xlabel(' f (Hz)'); ylabel('|H_2(f)|  (dB)');
maxy = max(20*log10(abs(H_2))) + 10; miny = maxy - 200;
axis([0 fs miny maxy]);
grid;
subplot(2,2,3);
plot(f_2048, angle(H_1))
xlabel('f (Hz)'); ylabel('Phase(H_1(f)) (rad)');
xlim([0 fs]);
grid;
subplot(2,2,4);
plot(f_2048, angle(H_2), 'Color', [0.8500, 0.3250, 0.0980])
xlabel(' f (Hz)'); ylabel('Phase(H_2(f)) (rad)');
xlim([0 fs]);
grid;

%As suggest in the assignment, we demodulate multiplying the modulated 
%signal by the previously extracted carriers
left_demod  = 2 .* y .* carry1 / A1;
right_demod = 2 .* y .* carry2 / A2;

%Cascade of lowpass and highpass as suggest in the assignment in order to 
%obtain the correct output
%Highpass notch; all the parameters are extract in a similar way we've
%done in the previous computation, knowing we're using a notch of the
%second order in [0, 20]
f_3db_notch = 20;
delta_notch = 2*pi*f_3db_notch/fs;
r_notch = 1 - delta_notch;
a1 = -2*r_notch;
a2 = r_notch*r_notch;

num_highpass = [1 -2 1];
den_highpass = [1 a1 a2];
H_1 = dfilt.df2t(num_highpass, den_highpass);
[H_highpass, w] = freqz(num_highpass, den_highpass, 'whole', 2048, fs);
H_highpass_dB = 20*log10(abs(H_highpass));

%Lowpass; in the 'ellip' specification we take the passband edge
%frequency of 8000Hz and set the other specifications as follow
[num_lowpass, den_lowpass] = ellip(15, 0.1, 100, 2*8000/fs, 'low');
H_2 = dfilt.df2t(num_lowpass, den_lowpass);
[H_lowpass, w] = freqz(num_lowpass, den_lowpass, 'whole', 2048, fs);
H_lowpass_dB = 20*log10(abs(H_lowpass));

%Plot of the highpass and lowpass filters, then the cascade
figure(4)                                  
subplot(2,1,1);
plot(f_2048, H_highpass_dB, 'Color', [0.4660, 0.6740, 0.1880])
title('Magnitude of the high pass filter');
xlabel(' f (Hz)'); ylabel('|H_{highpass}(f)| (dB)');
maxy = max(H_highpass_dB) + 10; 
miny = maxy - 200;
axis([0 fs miny maxy]);
grid;
subplot(2,1,2);
plot(f_2048, angle(H_highpass), 'Color', [0.4660, 0.6740, 0.1880])
title('Phase of the high pass filter');
xlabel(' f (Hz)'); ylabel('Phase(H_{highpass}(f)) (rad)');
xlim([0 fs]);
grid;

figure(5)                                   
subplot(2,1,1);
plot(f_2048, H_lowpass_dB, 'Color', [0.3010, 0.7450, 0.9330])
title('Magnitude of the low pass filter');
xlabel(' f (Hz)'); ylabel('|H_{lowpass}(f)|  (dB)');
maxy = max(H_lowpass_dB) + 10; 
miny = maxy - 200;
axis([0 fs miny maxy]);
grid;
subplot(2,1,2);
plot(f_2048, angle(H_lowpass), 'Color', [0.3010, 0.7450, 0.9330])
title('Phase of the low pass filter');
xlabel(' f (Hz)'); ylabel('Phase(H_{lowpass}(f)) (rad)');
xlim([0 fs]);
grid;

%Matlab implementation
H_cas = dfilt.cascade(H_1, H_2);
fvtool(H_cas);
title('Magnitude of the cascade of two filters')

y_dem = zeros(N, 2);
y_dem(:,1) = filter(H_cas, left_demod);
y_dem(:,2) = filter(H_cas, right_demod);

X_1 = fft(y_dem(:,1))/N;
X_2 = fft(y_dem(:,2))/N;
X_1_dB = 20*log10(abs(X_1));
X_2_dB = 20*log10(abs(X_2));

Y_dem = fft(y_dem);
Y_dem_norm = Y_dem/N;
Y_dem_dB = 20*log10(abs(Y_dem_norm));

%Manual implementation
% left_lowpass = filter(num_lowpass, den_lowpass, left_demod);
% left = filter(num_highpass, den_highpass, left_lowpass); 
% right_lowpass = filter(num_lowpass, den_lowpass, right_demod);
% right = filter(num_highpass, den_highpass, right_lowpass);
% y_dem = zeros(N, 2);
% y_dem(:,1) = left;
% y_dem(:,2) = right;

miny_1 = min(y_dem(:,1));
miny_2 = min(y_dem(:,2));
miny = max(miny_1, miny_2);
maxy_1 = max(y_dem(:,1));
maxy_2 = max(y_dem(:,2));
maxy = max(maxy_1, maxy_2);

maxY_1 = max(Y_dem_dB(:,1));
maxY_2 = max(Y_dem_dB(:,2));
maxY = max(maxY_1, maxY_2);

%Plot of x_1(nT) and x_2(nT)
figure(7)
subplot(2,1,1);
plot(f, X_1_dB, 'Color', [0.6350, 0.0780, 0.1840])
title('Spectrum of x_{1}');
xlabel(' f (Hz)'); ylabel('|X_1(f)| [dB]');
axis([0.1 d (maxY_1 - 250) maxY_1]);
xlim([0 fs]);
grid;
subplot(2,1,2);
plot(f, X_2_dB, 'Color', [0.6350, 0.0780, 0.1840])
title('Spectrum of x_{2}');
xlabel(' f (Hz)'); ylabel('|X_2(f)| [dB]');
axis([0.1 d (maxY_2 - 250) maxY_2]);
xlim([0 fs]);
grid;

%Plot of the new total spectrum in time domain
figure(8)              
plot(t, y_dem, 'Color', [0.6350, 0.0780, 0.1840])                                     
title('Demodulated signal: time domain');
xlabel('nT (s)'); ylabel('y(nT)');
axis([0.1 d miny maxy]);
grid;
%and in frequency domain
figure(9)
plot(f, Y_dem_dB, 'Color', [0.6350, 0.0780, 0.1840]) 
title('Demodulated signal: DFT'); 
xlabel('f (Hz)'); ylabel('|Y(f)| (dB)');
axis([0 fs (maxY - 250) maxY]);
grid;

%Comparison: before and after
figure(10)
grid on
plot(t, y, 'Color', [0.8500, 0.3250, 0.0980])
hold on 
plot(t, y_dem, 'Color', [0.9290, 0.6940, 0.1250])
title('Comparison between modulated signal and demodulated one in time');
legend({'Initial y(nT)','Final y(nT)'},'Location','northeast')
xlabel('nT (s)'); ylabel('y(nT)');
axis([0.1 d min(y) max(y)]);
hold off

audiowrite('Zabot_Aurora.wav', y_dem, fs);