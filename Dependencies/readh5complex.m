function [output] = readh5complex(filename)
%READH5COMPLEX reads a h5 file with \real and \imag datasets and returns a
%complex matrix
% Kylie Yeung 11/2024
if ~isempty(regexpi(filename,'\.h5$', 'once')), filename = filename(1:end-3); end

output_real=h5read(sprintf('%s.h5',filename), '/real');
output_imag=h5read(sprintf('%s.h5',filename), '/imag');

output=single(output_real+1i*output_imag);
end

