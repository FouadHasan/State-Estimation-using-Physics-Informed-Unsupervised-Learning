clear all
clc

load('StateEst_Case14.mat')
V=V(9001:10000,:)';
Z=Z(9001:10000,:)'/100;
clearvars Zri Vri Vma B G

n_scenario=1000;

%%
load('predicted_V.mat')
yr=real(Ybus);
yi=imag(Ybus);

for s=1:n_scenario
    s

    %Nodal injection measurement
    z=[real(Z(:,s));imag(Z(:,s))];
    
    %% Predicted Initialization
    vr=vr_pred(s,:)';
    vi=vi_pred(s,:)';
    v=vr+1j*vi;
    
    rmse_e=1;
    iter=0;
    
    while rmse_e>0.000001
        iter=iter+1;
        %% Getting nodal injection
        vr=real(v);
        vi=imag(v);
        
        p=((yr*vr-yi*vi).*vr+(yr*vi+yi*vr).*vi);
        q=((yr*vr-yi*vi).*vi-(yr*vi+yi*vr).*vr);
        S=[p;q]; %estimated nodal ijection
        
        % forming the jacobian
        dS_dvr=diag(v)*conj(Ybus)+conj(diag(Ybus*v));
        dS_dvi=1j*(conj(diag(Ybus*v))-diag(v)*conj(Ybus));
        
        H=[real(dS_dvr),real(dS_dvi);
            imag(dS_dvr),imag(dS_dvi)];
        
        A=H'*diag(ones(28,1))*H;
        b=H'*diag(ones(28,1))*(z-S);
        dVri=A\b;
        
        dVr=dVri(1:14,1);
        dVi=dVri(15:28,1);
        dV=dVr+1j*dVi;
        v_iter(:,iter)=v;
        v=v+dV;
        
        mse_e=mse(z,S);
        rmse_e=sqrt(mean((z-S).^2));
        
%         rmse_error(s,iter)=rmse_e;
        
        if iter ==500
            InfeasibilityFlag(s)=1;
            break;
        end        
    end%end of while
    iter
    itercount(s,1)=iter;
end%end of for

avg_iter=sum(itercount)/1000


% mse_e
% rmse_e
% [V v v*exp(-1j*angle(v(1)))]
% i

% for n=1:28
%     [N,edges] = histcounts(Vri(:,n),10);
%     maxIndex=find(N==max(N));
%     vri_ref(n,1)=(edges(maxIndex)+edges(maxIndex+1))/2;
% end