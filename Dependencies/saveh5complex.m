function saveh5complex(filename,input)
%SAVEH5COMPLEX saves complex input to a h5 file with \real and \imag datasets
% Kylie Yeung 11/2024
if ~isempty(regexpi(filename,'\.h5$', 'once')), filename = filename(1:end-3); end

h5create(sprintf('%s.h5',filename),'/real',size(input))
h5create(sprintf('%s.h5',filename),'/imag',size(input))

h5write(sprintf('%s.h5',filename),'/real',real(input))
h5write(sprintf('%s.h5',filename),'/imag',imag(input))
end

