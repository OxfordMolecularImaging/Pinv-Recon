function [F_z,k,dim] = getFz(k,mtx_acq, mtx_reco,npts)
%GETFZ get the reconstruction matrix in Z
    z = k(:,:,3);[~,i,~]=unique(z);
    Kz = single(z(sort(i)))./0.5*pi./mtx_reco(3).*mtx_acq(3);
    if mod(mtx_reco(3),2)==1
        ZM = single(-(mtx_reco(3)-1)/2:(mtx_reco(3)-1)/2); 
    else
        ZM = single(-(mtx_reco(3)/2):(mtx_reco(3)/2-1));
    end
    E_z = exp(1i*(-Kz(:)*ZM(:).'));
    F_z = pinv(E_z,0.5);

    clear z Kz ZM E_z
    
    % modify trajectory
    k = k(:,1:npts/mtx_acq(3),1:2);dim=2;
end

