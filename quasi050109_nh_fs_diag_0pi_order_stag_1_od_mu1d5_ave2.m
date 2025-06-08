%  N=987 jj=100 V=3.9 t=1.4e+02s
%  N=987 jj=100 V=4 t=1.5e+02s
% quasi022802_f_6_N_987_100realz_Delta_0d3

% 2025-3-2 15:57:03
% run('D:\OneDrive2\quasi_new_5_cluster\quasi022802_f_6_N_987_100realz_Delta_0d4.m')
% run('D:\OneDrive2\quasi_new_5_cluster\quasi022802_f_6_N_987_100realz_Delta_0d5.m')
%  N=987 jj=100 V=0 t=3.7s
% run('D:\OneDrive2\quasi_new_5_cluster\quasi022802_f_6_N_987_100realz_Delta_1d0.m')
%  N=987 jj=100 V=0.1 t=7.4s
%  set(gca,'YScale','log')
clear
f_str = 'f';
for kk=[1]
    tic
    
    format long
    J=1;
    BC=1;
    
    phi_list=2*pi*load('random_5000.dat');
    
    binary_f = load(['data1218_N10000_' f_str '.dat']);
    binary_f=binary_f';
    tic
    
    Deltas=[0.0];
    % Vs=[1.0];%dV=Vs(2)-Vs(1);
    Vs=[2.9:0.002:3.1];dV=Vs(2)-Vs(1);
    % 2025-3-25 23:14:43
    narray=[14:19]; %987
    % narray=[12:17]; %144:987
    % narray=[11:21]; %144:987 (11,89) (16,987)
    tic
    %  N=987 Delta=2 V=4 jj=1 t=105.3439
    %
    %  N=2584 jj=1 Delta=0.5 t=1.7s
    %  N=4181 jj=1 Delta=0.5 t=5.2s
    %  N=6765 jj=1 Delta=0.5 t=23s
    %
    %  N=6765 jj=1 Delta=0.1 t=24s
    %  N=10946 jj=1 Delta=0.1 t=1.8e+02s
    % quasi022802_f_3_xN_1pi_Delta_0d1_V_0d1
    % figure
    for iN=1:length(narray)
        
        Narray(iN)=fibonacci(narray(iN));
        N=Narray(iN);
        alpha=fibonacci(narray(iN)-1)/fibonacci(narray(iN));
        
        %for jj =20:100
        %     Deltas
        ave_ipr=zeros(length(Narray),length(Deltas));
        for iD=1:length(Deltas)
            Delta=Deltas(iD);
            for jj =1:100
                phi=phi_list(jj);
                for iV=1:length(Vs)
                    V=Vs(iV);
                    %disp([' N=' num2str(N) ' jj=' num2str(jj) ' Delta=' num2str(Delta) ' V=' num2str(V,1) ' t=' num2str(toc,2) 's'])
                    %disp([Delta V])
                    %alphas=V*cos(2*pi*alpha*(1:N)+phi)+Delta*(2*mod(1:N,2)-1); % 对称？
%                     h=0.1;
%                     %                     alphas=V*cos(2*pi*alpha*(1:N)+phi+1i*h); %+Delta*binary_f(1:N); % 对称？
%                     alphas=(V)*cos(2*pi*alpha*(1:N)+phi)+1i*h*(2*mod(1:N,2)-1); %+Delta*binary_f(1:N); % 对称？
%                     %H=-J*(diag(ones(1,N-1),1)+diag(ones(1,N-1),-1))+Deltas(iD)*diag(binary_f(1:N));
%                     %                 h=0.1;
%                     H=-J*(diag(ones(1,N-1),1)+diag(ones(1,N-1),-1))+diag(alphas);
%                     H(1,N)=-J*BC;
%                     H(N,1)=-J*BC;
                    
t=1;
     alphas=cos(2*pi*alpha*(1:N)+phi);
                    
                    mu=1.5;%Delta;
                    H=Vs(iV)*diag(alphas);
                    for j=1:N-1
                        H(j,j+1) = -t+mu*cos(2*pi*alpha*(j+0.5)+phi);
                        H(j+1,j) = -t+mu*cos(2*pi*alpha*(j+0.5)+phi);
                    end
                    H(1,N)= -t+mu*cos(2*pi*alpha*(N+0.5)+phi);
                    H(N,1)= -t+mu*cos(2*pi*alpha*(N+0.5)+phi);
            
            
                    [eigvec,eigval]=eig(H);
                    %psi=eigvec(:,1);
                    
                    % 对角化哈密顿量矩阵 H
                    [eigvec, eigval] = eig(H);
                    
                    % 提取本征值对角矩阵的对角元素（转换为列向量）
                    eigenvalues = diag(eigval);
                    
                    % 找到实部最小的本征值索引
                    [~, idx] = min(real(eigenvalues));
                    
                    % 选择对应的本征态
                    psi = eigvec(:, idx);
                    
                    %ipr(iV,iD,jj)=sum(abs(psi).^4);
                    ipr(iV,jj)=sum(abs(psi).^4);
                    if iV>1
                        f=abs(psi0'*psi);
                        %fs(iV,iN,iD,jj)=-2*log(f)/dV/dV;
                        fs(iV,jj)=-2*log(f)/dV/dV;
                    else
                        fs(iV,jj)=nan;
                    end
                    %                 fs_DN(iV,iD,iN,jj)=fs(iV,jj);
                    psi0=psi;
                end
                
                
                
                %             ave_ipr(iN,iD)=mean(ipr(iN,iD,:))
                disp([' N=' num2str(N) ' jj=' num2str(jj) ' Delta=' num2str(Delta) ' V=' num2str(V) ' t=' num2str(toc) 's'])
                
            end
            
            for iV=1:length(Vs)
                V=Vs(iV);
                data=fs(iV,:); % =fs(iV,jj);
                tfs(iV,iN) = prod(data.^(1/length(data)));
                %             tfs(iV) = prod(data)^(1/length(data));
            end
            
            figure
            
            for iN1=1:iN %length(narray)
                
                Narray(iN1)=fibonacci(narray(iN1));
                plot(Vs,tfs(:,iN1),'.')
                hold on
            end
            %
            %             figure
            %         plot(Vs,ave_fs1,'o')
            xlabel('V');
            ylabel('\chi_F');
            set(gca,'FontSize',15);
            title([ 'stag N=' num2str(N)  ' V=' num2str(V) ' mu=' num2str(mu) ' phi=' num2str(phi)])
            disp_legend_Narray% hold on
            currentFileName = mfilename;
            %             save([currentFileName '.mat']);
            
% 将 y 轴改为对数尺度
set(gca, 'YScale', 'log');  % gca 表示获取当前坐标轴
    hold on;  % 保持图形叠加
    
    % 绘制竖直线
    x_val = 2*exp(-0.1);
%     xline(x_val, '--k', 'LineWidth', 0.5);
    % fs(1,:)=nan;
%             ylim([0 50])
%             saveas(gcf, [currentFileName  '_fs_N_' num2str(N) ' .png']);
            saveas(gcf, [currentFileName  '_fs.png']);
    % 将数组保存为文本格式的 .dat 文件（空格分隔）
% save([currentFileName  '_fs_N_' num2str(N) '.dat'], 'tfs', '-ascii', '-double');
save([currentFileName  '_fs.dat'], 'tfs', '-ascii', '-double');
        end
        % %     hold on
        
        %         saveas(gcf, [currentFileName  '_ave.png']);
        
        %         save([currentFileName '_ipr_N_' num2str(N) '.dat']', 'fs', '-ascii');
        
        %     close all
        % fs_DN(iV,iD,iN,jj)
        % figure
        %     ave_fs=mean(fs_DN,4); %(iV,iD,iN,jj);
        %     % plot(Vs, ave_fs,'o')
        %     for iN=1:length(narray)
        %         for iV=1:length(Vs)
        %             for iD=1:length(Deltas)
        %                 data=fs_DN(iV,iD,iN,:); % =fs(iV,jj);
        %                 gmean_fs(iV,iD,iN) = prod(data)^(1/length(data));            % 计算乘积并取 n 次方根
        %             end
        %         end
        %     end
        
        %
        %     currentFileName = [mfilename '_' f_str '_' '_Nmax_' num2str(N)];
        %     %     %save([currentFileName '.mat']);
        %     %     % end
        %     %     iD=1;
        %     %     ave_fs1=reshape(ave_fs(:,iD,:),length(Vs),length(Narray));
        %     %     gmean_fs1=reshape(gmean_fs(:,iD,:),length(Vs),length(Narray));
        %     %     currentFileName = [mfilename '_N_' num2str(N)];
        %
        %     save([currentFileName '_Vs.dat']', 'Vs', '-ascii');
        %     %save([currentFileName '_ave_fs1.dat']', 'ave_fs1', '-ascii');
        %     save([currentFileName '_tfs.dat']', 'tfs', '-ascii');
    end
    figure
    
    for iN=1:length(narray)
        
        Narray(iN)=fibonacci(narray(iN));
        plot(Vs,tfs(:,iN),'.')
        hold on
    end
    %
    %             figure
    %         plot(Vs,ave_fs1,'o')
    xlabel('V');
    ylabel('\chi_F');
    set(gca,'FontSize',15);
    title([ ' N=' num2str(N)  ' V=' num2str(V) ' g=' num2str(h) ' phi=' num2str(phi)])
    disp_legend_Narray% hold on
    hold on;  % 保持图形叠加
    
    % 绘制竖直线
    x_val = 2*exp(-0.1);
    xline(x_val, '--r', 'LineWidth', 1.5);
    % fs(1,:)=nan;
    %
    % figure
    % plot(Narray, ipr,'o');
    % ylabel('P');
    % xlabel('N');
    % set(gca,'FontSize',15);
    % title([ ' V=' num2str(V) ' Delta=' num2str(Delta) ' phi=' num2str(phi)])
    %
    %
    % figure
    % plot(log(Narray), log(ipr),'o');
    % ylabel('log P');
    % xlabel('log N');
    % set(gca,'FontSize',15);
    % title([ ' V=' num2str(V) ' Delta=' num2str(Delta) ' phi=' num2str(phi)])
    %
    % end
    
    for x=[]
        % % save([currentFileName '_err_ipr.dat']', 'err_ipr', '-ascii');
        figure
        plot(Vs,ave_fs1,'o')
        xlabel('V');
        ylabel('\chi_F');
        set(gca,'FontSize',15);
        title([ ' N=' num2str(N)  ' V=' num2str(V) ' Delta=' num2str(Deltas(iD)) ' phi=' num2str(phi)])
        disp_legend_Narray% hold on
        saveas(gcf, [currentFileName  '_ave.png']);
        figure
        % plot(Vs,ave_fs1,'o')
        plot(Vs,gmean_fs1,'o')
        xlabel('V');
        % ylabel('\chi_F');
        ylabel('\chi_{F,typ.}');
        set(gca,'FontSize',15);
        title([ ' N=' num2str(N)  ' V=' num2str(V) ' Delta=' num2str(Deltas(iD)) ' phi=' num2str(phi)])
        disp_legend_Narray% hold on
        saveas(gcf, [currentFileName '_typ.png']);
        % figure
        % plot(Vs,gmean_fs,'o')
        % ylabel('\chi_{F,typ.}');
        % xlabel('V');
        % set(gca,'FontSize',15);
        % title([ ' N=' num2str(N)  ' V=' num2str(V) ' Delta=' num2str(Delta) ' phi=' num2str(phi)])
        % disp_legend_Narray% hold on
        % saveas(gcf, [currentFileName '_typfs.png']);
        
        % ave_fs=mean(fs,2);
        % plot(Vs, ave_fs,'o')
        % plot(Vs,gmean_fs,'o')
        % plot(Vs, fs,'o');
        % save([currentFileName '_Vs.dat']', 'Vs', '-ascii');
        % save([currentFileName '_.dat']', '', '-ascii');
        % % ave_ipr=0
        % for iV=1:length(Vs)
        %     errorbar(1./(Narray), (ave_ipr(:,iV)),err_ipr(:,iV),'o');
        %     hold on
        % end
        % ave_ipr=mean(ipr,2);
        %
        % for iN=1:length(narray)
        %     %     N=Narray(iN);
        %     for iV=1:length(Vs)
        %         %         V=Vs(iV);
        %         err_ipr(iN,iV) = std(ipr(iN,iV,:))/sqrt(N);
        %     end
        % end
        %
        % % plot(1./(Narray(1:7)), (iprave(1:7,:)),'o');
        % % plot(1./(Narray(:)), (ave_ipr(:,:)),'o');
        % hold on
        % plot(Vs, (ave_ipr(:,:)),'x');
        % ylim([0 0.1])
        % figure
        % plot(log10(Narray(:)), log10(ave_ipr(:,:)),'o');
        % ylabel('log P');
        % xlabel('log N');
        % set(gca,'FontSize',15);
        % title([ ' V=' num2str(V) ' Delta=' num2str(Delta) ' phi=' num2str(phi)])
        % ylim([0 1])
        
        %
        % set(gca, 'YScale', 'log'); % 设置当前坐标轴的 Y 轴为对数刻度
        % grid on;                    % 添加网格
        
        % figure
        % plot(1./(narray), (ipr),'o');
        % ylabel('P');
        % xlabel('1/n');
        % set(gca,'FontSize',15);
        % title([ ' V=' num2str(V) ' Delta=' num2str(Delta) ' phi=' num2str(phi)])
        
        
        % disp_legend_Narray
        
        % figure
        % ps=20;
        % ipr_ave=mean(ipr,3);
        % [X,Y]=meshgrid(Vs,Deltas);Z=real(ipr_ave');
        % %  [X,Y]=meshgrid(aarray,harray);yarray=real(gap);
        % %yarray=real(d2E0_trans');
        % meshplotvar=pcolor(X,Y,Z);
        % ps=15;
        % % ylabel('h','fontsize',ps);xlabel('\alpha','fontsize',ps);
        % xlabel('V','fontsize',ps);
        % % ylabel('\alpha','fontsize',ps);
        % ylabel('\Delta','fontsize',ps);
        % %
        % box on
        % set(meshplotvar,'edgecolor','none','facecolor','interp');
        %
        % set(get(gca,'XLabel'),'Fontsize',ps)
        % set(get(gca,'YLabel'),'Fontsize',ps)
        % set(gca,'FontSize',ps);
        % box on
        % colorbar
        % % title('N=987;alpha=610/987;')
        % title([ ' N=' num2str(N)  ' phi.ave=' num2str(phi)])
        % 获取当前文件的文件名（不带路径和扩展名）
        
        %
        % currentFileName = mfilename;
        % % title(currentFileName);
        % disp(currentFileName);
        % % %
        % % save([currentFileName '.dat']', 'mipr', '-ascii');
        % % save([currentFileName '_ipr.dat']', 'ipr', '-ascii');
        % % save([currentFileName '_ipr.dat']', 'ipr', '-ascii');
        % % save([currentFileName '_ave_ipr.dat']', 'ave_ipr', '-ascii');
        % % save([currentFileName '_err_ipr.dat']', 'err_ipr', '-ascii');
        % save([currentFileName '_Narray.dat']', 'Narray', '-ascii');
        % save([currentFileName '_Deltas.dat']', 'Deltas', '-ascii');
        % % saveas(gcf, 'quasi122701_nh_3d_4.png');
    end
end
