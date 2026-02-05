% ============================================================
% Reproduce Fig.1(a): IPR vs correlation strength c
% PRB 106, L220201 (2022)
% Vallejo-Fabila & Torres-Herrera
% ============================================================

clear;% clc;

%% ------------------ Parameters ------------------
L = 14;                  % system size
hList = [0.5, 3.75, 6];  % disorder strengths
cList = logspace(-1,1.3,40);   % c from ~0.1 to ~20
numEig = 100;            % eigenstates closest to E=0
runAvgWin = 10;          % running average window

rng(1);                  % fixed seed for reproducibility

%% ------------------ Build basis Sz = 0 ------------------
fprintf('Building Sz=0 basis...\n');
basis = [];
for s = 0:(2^L-1)
    if sum(bitget(s,1:L)) == L/2
        basis(end+1) = s; %#ok<SAGROW>
    end
end
D = length(basis);
fprintf('Hilbert space dimension D = %d\n',D);

stateIndex = containers.Map(basis,1:D);

%% ------------------ Spin operators ------------------
Sx = [0 1;1 0]/2;
Sy = [0 -1i;1i 0]/2;
Sz = [1 0;0 -1]/2;

%% ------------------ Loop over h ------------------
figure; hold on;
colors = {'r','g','b'};
tic
for hh = 1:length(hList)
    h = hList(hh);
    IPRc = zeros(length(cList),1);
disp([' h=' num2str(h) ' t=' num2str(toc) ])
    for ic = 1:length(cList)
        c = cList(ic);
disp([' h=' num2str(h) ' c=' num2str(c) ' t=' num2str(toc) ])

        % ---- generate correlated disorder hk ----
        hk = correlated_disorder(L,h,c);

        % ---- build Hamiltonian in Sz=0 sector ----
        H = zeros(D,D);

        for ii = 1:D
            s = basis(ii);

            for site = 1:L
                j = mod(site,L)+1;

                si = bitget(s,site);
                sj = bitget(s,j);

                % Sz Sz
                zi = (si==1)*0.5 + (si==0)*(-0.5);
                zj = (sj==1)*0.5 + (sj==0)*(-0.5);
                H(ii,ii) = H(ii,ii) + zi*zj;

                % Sx Sx + Sy Sy
                if si ~= sj
                    sf = bitset(bitset(s,site,sj),j,si);
                    jj = stateIndex(sf);
                    H(ii,jj) = H(ii,jj) + 0.5;
                end
            end

            % disorder term
            for site = 1:L
                si = bitget(s,site);
                zi = (si==1)*0.5 + (si==0)*(-0.5);
                H(ii,ii) = H(ii,ii) + hk(site)*zi;
            end
        end

        % ---- diagonalize ----
        [V,E] = eig(H);
        E = diag(E);

        % ---- select eigenstates near E=0 ----
        [~,idx] = sort(abs(E));
        idx = idx(1:numEig);

        % ---- compute IPR ----
        ipr = 0;
        for k = idx'
            psi = V(:,k);
            ipr = ipr + sum(abs(psi).^4);
        end
        IPRc(ic) = ipr/numEig;
    end

    % ---- running average ----
    IPRsmooth = movmean(IPRc,runAvgWin);

    loglog(cList,IPRsmooth,colors{hh},'LineWidth',2);
end

% ---- reference lines ----
yline(3/(D+2),'k--','GOE');
yline(1,'c--','Poisson');

xlabel('c');
ylabel('IPR');
legend('h=0.5','h=3.75','h=6','Location','southwest');
set(gca,'FontSize',12);
box on;

title('Reproduction of Fig.1(a): IPR vs correlation strength c');

%% ================== FUNCTIONS ==================
function hk = correlated_disorder(L,h,c)
    V = rand(1,L+1);
    X = zeros(1,L);

    X(1) = V(1) + V(2)/c;
    for k = 2:L
        X(k) = cdf_sum_uniform(X(k-1),c) + V(k+1)/c;
    end

    F = cdf_sum_uniform(X,c);
    hk = h*(2*F - 1);
end

function F = cdf_sum_uniform(x,c)
    F = zeros(size(x));
    for i = 1:length(x)
        xi = x(i);
        if xi < 0
            F(i) = 0;
        elseif xi < 1/c
            F(i) = c*xi^2/2;
        elseif xi < 1
            F(i) = xi - 1/(2*c);
        elseif xi < 1 + 1/c
            F(i) = 1 - c*(1 + 1/c - xi)^2/2;
        else
            F(i) = 1;
        end
    end
end
