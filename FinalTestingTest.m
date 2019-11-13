clear all;
clc;
offset = 4;
Gd = 3; Gr = 2; Td = 1; Tr = 2;
rows = 30;
cols = 25;
x = rand(rows,cols);
x(5,11)= 100;
x(5,12)= 121;
x(5,13)= 124;
sig_fft2 = fft2(x,rows,cols);
sig_fft2 = sig_fft2(1:rows,1:cols);    % given rows/2 in the code.. this make array out of bound error.. look into it once again. 
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

GrdGridArray = [];
threshold_cfar = [];  % Threshold values Vector. 
signal_cfar = [];   % Final Signal Vector
for i = Gr+Tr+1:rows-Gr-Tr
for j = Gd+Td+1 : cols-Gd-Td
    CUT = RDM(i,j);
    FullGrid = RDM(i-Gr-Tr: i+Gr+Tr, j-Gd-Td: j+Gd+Td);
    GrdGrid = RDM(i-Gr : i+Gr,j-Gd : j+Gd);
    GrdGridArray = [GrdGridArray,{GrdGrid}];
    TrailBitSum = sum(sum(db2pow(FullGrid))') - sum(sum(db2pow(GrdGrid))'); 
    % Can also be written as 
% TrailBitSum = sum(FullGrid,'all') - sum(GrdGrid, 'all'); 
    TrailBits = numel(FullGrid)- numel(GrdGrid);
% or TrailBits1 = prod(size(FullGrid)) - prod(size(GrdGrid));
% or (2Gr+2Tr+1)(2Gd+2Td+1)-(2Gr+1)(2Gd+1)
    threshold = TrailBitSum/TrailBits;
    threshold_db = pow2db(threshold);
    threshold_scaled = threshold_db + offset;
    threshold_cfar = [threshold_cfar, {threshold_scaled}];
    signal_level = CUT-threshold_scaled;
    if (signal_level<0)
          CUT = 0;
    else
        CUT = 1;
    end
%         signal_cfar = [signal_cfar,{CUT}];
    

end
end

doppler_axis = linspace(-100,100,cols);
range_axis = linspace(-200,200,rows/2)*((rows/2)/400);
figure,surf(doppler_axis,range_axis,signal_cfar);
colorbar;
