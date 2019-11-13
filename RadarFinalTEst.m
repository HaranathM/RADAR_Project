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

init_range = 40;         % initial position(relative to radar) in meters. 
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
x1 = t(2)-t(1);
for i=1:length(t)         
    
    %x = Nd*Tchirp/((Nr*Nd)-1)    This is equivalent to t(2)-t(1)      
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = init_range + (init_vel * t(i));      % Refer this again!!!
    % delayed time or trip time
  td(i) = 2*r_t(i)/c;
%     td = 2*init_range/c;
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i)  =  cos(2*pi*(fc*i+(slope*(t(i)^2)/2)));
    Rx(i)  =  cos(2*pi*(fc*(t(i)-td(i))+(slope*((t(i)-td(i))^2)/2)));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix1(i) = Tx(i).*Rx(i);
    
end
%     Mix1 = Tx.*Rx;

%% RANGE MEASUREMENT

% figure(1)
% plot(1000*t(1:50) ,Tx(1:50));
% figure(2)
% plot(5000*t(1:100) ,Rx(1:100));
 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
% Mix_mat = vec2mat(Mix,Nd)   ---    This needs Communication Toolbox though. 
Mix_mat = reshape(Mix1,[Nr,Nd]);
% figure(3);
% plot(1000*t(1:50) ,Mix1(1:50));
 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.

% fft_mix_trans = fft2(Mix_mat');
% fft_mix = fft_mix_trans';

fft_mix1 = fft(Mix_mat,[],1)/Nr;

% *%TODO* :
% Take the absolute value of FFT output
fft_mix = abs(fft_mix1);
 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
 fft_mix_half = fft_mix(1:Nr/2+1);


%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)


 % *%TODO* :
 % plot FFT output 
 fs = 1/x1;
f = (/length(fft_mix_half))*(0:(1024/2));
%  /1028*128
 plot(f,fft_mix_half);

 
% axis ([0 500 0 0.3]);

