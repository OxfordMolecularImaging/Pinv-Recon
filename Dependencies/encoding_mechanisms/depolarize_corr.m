function [d_corr] = depolarize_corr(d,flipvec)
%DEPOLARIZE_CORR Summary of this function goes here
%   Detailed explanation goes here
if ~exist("flipvec","var")
    load("C:\Users\kylie\OneDrive - Nexus365\LEARN\Waveforms\xenon_vent\Xenon_dualres_excs766_v3_vap_vals.mat")
    warning("Flip vec not defined. Assume to be as Xenon_dualres_excs766_v2_vap_vals.mat")
else
    Flip_vent=flipvec;
end
cos_factor=cumprod(cosd(Flip_vent));
b1_factor=repmat([1 cos_factor(1:end-1)]',1,size(d,2));
d_corr=d./b1_factor;
end

