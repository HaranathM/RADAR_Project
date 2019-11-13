% HAve to calculate velocities of the targets with the respective
% doppler frequency shifts for the given data
fd = [3e3 4.5e3 11e3 -3e3]     % Given doppler frequency shifts in Hertz
ft = 77e9;  % Radar operating frequency in Hz (77GHz)
c = 3e8;    % Speed of light

% Overview of the solution : fd = 2*vr/lambda; lambda = c/ft

lambda = c/ft;
rel_velocity = fd*lambda/2     % velocity in m/s
%disp(rel_velocity);

% positive velocities corresponds to approaching targets while 
% negative velocities corresponds to receding targets