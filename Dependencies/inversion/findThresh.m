function [F,pinv_thresh,pinv_time] = findThresh(E,svdE,svd_threshold,plt_svd)
% define pinv threshold relative to svd vals
cumulative_energy = cumsum(svdE) / sum(svdE);
g = find(cumulative_energy >= svd_threshold, 1);
if isempty(g);g=size(svdE,1);end
pinv_thresh=svdE(g);
if plt_svd; figure;plot(svdE);hold on;set(gca,'YTick',gather(pinv_thresh));grid on; end

pinv_start=tic;
F=pinv(E,pinv_thresh);
pinv_time=toc(pinv_start);
end

