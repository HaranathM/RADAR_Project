%% Defining a sinusoidal signal with the following parameters
Fs = 1000;    %Sampling Freq
T = 1/Fs;     % Sample width
L = 1500;     % Length of the signal
t = (0:L-1)*T;  % Time vector

%% Form a signal containing a 77 Hz sinusoid of amplitude 0.7 and a 
% 43Hz sinusoid of amplitude 2.

signal = 0.7* cos (2*pi*77*t)+ 2* cos(2*pi*43*t);

% Adding noise to the signal
n_signal = signal + 2*randn(size(t));

% Let us plot the noisy signal in the time domain. 
% We can observe that it's difficult to identify the  
% frequency components by looking at the signal n_signal(t). 

% figure(1)
% plot(t(1:50),n_signal(1:50))
% figure(2)

% Both the representations are same except that the time is represented in 
% seconds in the below case, rather than a very complicated decimal number
% as in the above case. Decomment above 3 lines of code to get clarity. 

plot(1000*t(1:50) ,n_signal(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('n_signal')

%% Compute the Fourier transform of the signal. 
freq_rep_signal = fft(n_signal)  

% Compute the two-sided spectrum P2. 
% Then compute the single-sided spectrum P1 
% based on P2 and the even-valued signal length L.

P2 = abs(freq_rep_signal);
P1 = P2(1:L/2);

%% Plotting

% Let's normalize the frequencies before plotting the fft.  
f = Fs*(0:(L/2)-1)/L;  

figure(3)
plot(f,P1); 
title('Single-Sided Amplitude Spectrum of X(t)');
xlabel('f (Hz)');
ylabel('|P1(f)|');








