function [cMap] = RGcMap()
%RGCMAP gives a colormap that is green for larger than 1 and red for
%smaller than 1
    greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
    redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
    cMap = [redColorMap; greenColorMap; zeros(1, 256)]';
end

