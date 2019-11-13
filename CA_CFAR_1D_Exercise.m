% Steps to implement CFAR Algorithm: 
% 1. Define the number of training cells and guard cells
% 2. Start sliding the window one cell at a time across the complete FFT 1D array. Total window size should be: 2(T+G)+CUT
% 3. For each step, sum the signal (noise) within all the leading or lagging training cells
% 4. Average the sum to determine the noise threshold
% 5. Using an appropriate offset value scale the threshold
% 6. Now, measure the signal in the CUT, which is T+G+1 from the window starting point
% 7. Compare the signal measured in 5 against the threshold measured in 4
% 8. If the level of signal measured in CUT is smaller than the threshold measured, then assign 0 value to the signal within CUT.

close all;
% Number of data points
Ns = 1000;
s = randn(Ns,1);     %Random noise
s = abs(s)     % See the difference when absolute value is not taken.(Its drastic!!) 
s([100 200 300 700]) = [8 9 4 11];   % Let the targets be at these points 
% 99, 200, 300, 700 with respective amplitudes 8 9 4 11.
% plot(s);
% fft(s);
% figure(2)
% plot(s) 

% The plot above is same as the plot above What I mean is the time 
% and frequency domain responses have no distinction in this case
% as the signal is random signal which is not varying with time in
% the first case. 
% In other words, there's no need to calculate the FFT for the given 
% signal(as in step 2) as the signal itself represents the one in frequency 
% domain. 

% Lets apply the above steps
% Define Training and Guard Cells

G = 4;   %  Guard Cell count
T = 10;   %  Training Cell count
offset = 4; % Adds room above noise threshold for desired SNR
threshold_cfar = [];  % Threshold values Vector. 
signal_cfar = [];   % Final Signal Vector

for i = 1: Ns-2*G-2*T 
    CUT = s(i+G+T);   % Array indexing starts from 1 (not 0)
    signal_sum = sum(s(i+G+T+1:i+2*G+2*T)); % Sum of noise in all leading training cells

% for i = 1: Ns-G-T 
%     CUT = s(i+G+T);   % Array indexing starts from 1 (not 0)
%     signal_sum = sum(s(i:i+T-1));  % Sum of noise in all lagging training cells 
    threshold = signal_sum/T;
    threshold_scaled = threshold*offset;
    threshold_cfar = [threshold_cfar, {threshold_scaled}];
    signal_level = CUT-threshold_scaled;
    if (signal_level<0)
          CUT = 0;        
    end
        signal_cfar = [signal_cfar,{CUT}];
         
end

% figure(2);
plot(cell2mat(signal_cfar),'g--');    

figure,plot(s);
hold on,plot(cell2mat(circshift(threshold_cfar,G)),'r--','LineWidth',2)
hold on, plot (cell2mat(circshift(signal_cfar,(T+G))),'g--','LineWidth',4);
legend('Signal','CFAR Threshold','detection')
    
    
    
    
    
    
