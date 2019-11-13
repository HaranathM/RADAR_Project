% Observation : As beat frequency increases, the Range will increase
% i.e., the target is moving far away from the Radar.

% Defining the needed parameters

fb = [0 1.1 13 24]; %beat frequencies in MHz
%Can also be defined as fb = [0 1.1e6 13e6 24e6]
c = 3*10^8;  % velocity of light m/s
dres = 1;   %Range Resolution in meters
Rmax = 300; %Value in meters

%Calculating Chirp BW and Sweep Time
% Sweep time is generally 5 to 6 times the time taken for a round trip.
Bsweep = c/(2*dres);
Tsweep = 5.5*2*Rmax/c;

% Range Calculation in meters

Range = c*Tsweep*fb*10^6/(2* Bsweep)
%disp(Range)
