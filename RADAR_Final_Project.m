clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
% Velocity Resolution =  3m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant

init_range = 160;         % initial position(relative to radar) in meters. 
init_vel = 50;            % initial velocity in m/s.
 


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

c = 3e8;       % speed of light (m/s)
Rmax = 200;    % Max range in meters.
dres = 1;      % Range resolution in meters.     
B = c/(2*dres); 
Tchirp = 5.5*2*Rmax/c;
slope = B/Tchirp;

%Operating carrier frequency of Radar 
fc= 77e+9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 
x1 = t(2)-t(1)
for i=1:length(t)         
    
    %x = Nd*Tchirp/((Nr*Nd)-1)    This is equivalent to t(2)-t(1)      
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = init_range + init_vel * t(i);      % Refer this again!!!
    % delayed time or trip time
    td(i) = 2*r_t(i)/c;
    % *%TODO* :
    %For each time sample we need to update the transmitted and
    %received signal. 
    Tx(i)  =  cos(2*pi*(fc*i+(slope*(t(i)^2))/2));
    Rx(i)  =  cos(2*pi*(fc*(t(i)-td(i))+(slope*((t(i)-td(i))^2))/2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix1(i) = Tx(i).*Rx(i);
    
end
%     Mix1 = Tx.*Rx;

%% RANGE MEASUREMENT

%% Try along these lines - Begin section (Aditya)
Mix_mat = reshape(Mix1,[1, Nr*Nd]);
fft_mix = fft(Mix_mat, Nr, 2);
P2 = abs(fft_mix/(Nr*Nd));
P1 = P2(1:Nr/2+1);

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

% Plotting
f = fc*(0:(Nr/2))/(Nr);
plot(f, P1);
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
%% - End Section (Aditya)


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
% Mix_mat = vec2mat(Mix,Nd)   ---    This needs Communication Toolbox though. 
Mix_mat = reshape(Mix1,Nr,Nd);
 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.

% fft_mix_trans = fft2(Mix_mat');
% fft_mix = fft_mix_trans';

fft_mix = fft(Mix_mat,[],2);

% *%TODO* :
% Take the absolute value of FFT output
fft_mix = abs(fft_mix);
 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
 fft_mix_half = fft_mix(1:1024*128/2);


%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)


 % *%TODO* :
 % plot FFT output 
 
 f = fc*(0:(1024*128/2)-1)/1024*128;
%  plot(f,fft_mix_half);
 
 plot (fft_mix_half)

 
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix1,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 2;
Td = 4;

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 1;
Gd = 2;
CUT = RDM(Gr+Tr+1,Gd+Td+1)
% *%TODO* :
% offset the threshold by SNR value in dB
offset_dB = 6;

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR
   
   GrdGridArray = [];
threshold_cfar = [];  % Threshold values Vector. 
signal_cfar = [];   % Final Signal Vector
for i = Gr+Tr+1:Nr/2-Gr-Tr
for j = Gd+Td+1 : Nd-Gd-Td
    CUT = RDM(i,j);
    FullGrid = RDM(i-Gr-Tr: i+Gr+Tr, j-Gd-Td: j+Gd+Td);
    GrdGrid = RDM(i-Gr : i+Gr,j-Gd : j+Gd);
    GrdGridArray = [GrdGridArray,{GrdGrid}];
    TrailBitSum = sum(sum(db2pow(FullGrid))') - sum(sum(db2pow(GrdGrid))'); 
    % Can also be written as 
% TrailBitSum = sum(db2pow(FullGrid),'all') - sum(dB2pow(GrdGrid), 'all'); 
    TrailBits = numel(FullGrid)- numel(GrdGrid);
% or TrailBits1 = prod(size(FullGrid)) - prod(size(GrdGrid));
% or (2Gr+2Tr+1)(2Gd+2Td+1)-(2Gr+1)(2Gd+1)
    threshold = TrailBitSum/TrailBits;
    threshold_dB = pow2db(threshold);
    threshold_scaled = threshold_dB + offset_dB;    
    threshold_cfar = [threshold_cfar, {threshold_scaled}];
    signal_level = CUT-threshold_scaled;
    if (signal_level>0)
          CUT = 1;
    else CUT = 0;
    end
        signal_cfar = [signal_cfar,CUT];
    

end
end
   

signal_cfar=reshape(signal_cfar,[Nr/2,Nd]);



% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 








% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);

figure,surf(doppler_axis,range_axis,signal_cfar);
colorbar;


 
 
