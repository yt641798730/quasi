% xtt !!!
for xtt=[]
    clear
    tic
    
    format long
    J=1;
    BC=1;
    
    % load phi
    phi_list=2*pi*load('random_5000.dat');
    binary_f = load('data1218_N10000_f.dat');
    binary_f=binary_f';
    tic
    
    Deltas=[0:0.04:2];
    Vs=[0:0.08:4];dV=Vs(2)-Vs(1);
    % narray=[16]; %987
    narray=[16]; %144:987
    tic
    %  N=987 Delta=2 V=4 jj=1 t=105.3439
    %  N=987 Delta=-1 V=-1 jj=39 t=3.1e+04s
    %  N=987 Delta=-1 V=-1 jj=40 t=3.2e+04s
    %  N=987 Delta=-1 V=-1 jj=41 t=3.2e+04s
    %  N=987 Delta=-1 V=-1 jj=42 t=3.3e+04s
    %  N=987 Delta=-1 V=-1 jj=43 t=3.4e+04s
    for iN=1:length(narray)
        
        Narray(iN)=fibonacci(narray(iN));
        N=Narray(iN);
        alpha=fibonacci(narray(iN)-1)/fibonacci(narray(iN));
        
        %for jj =20:100
        for jj = 1:100 %% perfect !!!
            phi=phi_list(jj);
            Delta=-1;V=-1;
            disp([' N=' num2str(N) ' Delta=' num2str(Delta) ' V=' num2str(V) ' jj=' num2str(jj) ' t=' num2str(toc,3) 's'])
            for iD=1:length(Deltas)
                Delta=Deltas(iD);
                
                for iV=1:length(Vs)
                    V=Vs(iV);
                    %alphas=V*cos(2*pi*alpha*(1:N)+phi)+Delta*(2*mod(1:N,2)-1); % 对称？%% perfect !!!
                    alphas=V*cos(2*pi*alpha*(1:N)+phi)+Delta*binary_f(1:N); % 对称？ %% perfect !!!
                    H=-J*(diag(ones(1,N-1),1)+diag(ones(1,N-1),-1))+diag(alphas);
                    H(1,N)=-J*BC;
                    H(N,1)=-J*BC;
                    [eigvec,eigval]=eig(H);
                    %  for kk=1:length(N) % psi=eigvec(:,kk);ipr(iV,iD,kk)=sum(abs(psi).^4);   end
                    psi=eigvec(:,1);
                    ipr(iV,iD)=sum(abs(psi).^4 );
                    if iV>1
                        f(iV,iD)=abs(psi0'*psi);%% perfect !!!
                        %fs(iV,iN)=-2*log(f)/dV/dV;
                    else
                        f(iV,iD)=1; %fs(iV,iN)=nan;
                    end
                    if iV==1
                        psi0=psi;
                    end
                    
                    fj(iV,iD,jj)=f(iV,iD);
                    iprj(iV,iD,jj)=ipr(iV,iD);
                end
            end
            
            
            for iD=1:length(Deltas) %% perfect !!!
                for iV=1:length(Vs)
                    %fj(iV,iD)=mean(fj(iV,iD,:));
                    fjj=fj(iV,iD,1:jj);
                    ave_f(iV,iD) = mean(fjj);
                    err_f(iV,iD) = std(fjj)/sqrt(N);
                    
                    iprjj=iprj(iV,iD,1:jj);
                    ave_ipr(iV,iD) = mean(iprjj);
                    err_ipr(iV,iD) = std(iprjj)/sqrt(N);
                end
            end
            currentFileName = mfilename;%% perfect !!!%% perfect !!!
            save([currentFileName '_tmp.mat']);
        end
        
        
    end
end
load('quasi030401_f_9_fs_phase_N_987_phi_many_V_1_ipr_redo_stag_tmp.mat')
currentFileName = mfilename;
save([currentFileName '.mat']);
figure
ps=20;
[X,Y]=meshgrid(Vs,Deltas);Z=real(ave_ipr'); %gap);d2E0_trans');
meshplotvar=pcolor(X,Y,Z);%% perfect !!!
ps=15;
xlabel('V','fontsize',ps);
ylabel('\Delta','fontsize',ps);
box on
set(meshplotvar,'edgecolor','none','facecolor','interp');
set(get(gca,'XLabel'),'Fontsize',ps)
set(get(gca,'YLabel'),'Fontsize',ps)
set(gca,'FontSize',ps);
box on
colorbar
title([ 'ipr N=' num2str(N)  ' phi=' num2str(phi)])

currentFileName = mfilename;
disp(currentFileName);
save([currentFileName '_ave_f.dat']', 'ave_f', '-ascii');
save([currentFileName '_ave_ipr.dat']', 'ave_ipr', '-ascii');
saveas(gcf, [currentFileName '_ave_ipr.png']);

