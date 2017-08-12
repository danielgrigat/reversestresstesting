function [R1] = reverse_policy_SREP(source,xT1,pol,e1,p,varargin)
%% equity capital policy plotting

%input
%source     = for STOXX variable density (<1); (always year 2015)
%               for full STOXX: (1);
%xT1        = desired target state
%pol        = policy: 1 relative, 2 K_i based
%e1         = increase of equity for K_i based policy (pol=2); e.g. for 5% enter 0.05
%p          = plot on(1) off(0)
%varargin{1}=   when source <1 is iterator for density simulation from 1 to 100
%               when source =1 is indicator of year


%% set up

xT=0.1:0.1:3; %target state
adju=0.1:0.1:3; %\lambda_max
T=20;
ind0=find(xT==xT1); %find location of result for specific target state

%% load necessary data

if source < 1 %STOXX changing density

    [~, ~, ~, ~, ~, equity] = import_stoxx('2015');
    n=length(equity);
    
    so = int64(source*10);
    load(['2015/density/adj_T20_d' num2str(so) '.mat']);
    adj=adj(:,:,varargin{1});
    load(['2015/density/u.mat']); 
    u_T20=u(:,:,:,varargin{1},so);   
    load(['2015/density/E_node.mat']); 
    E_node_T20=E_node(:,:,varargin{1},so);   
    
elseif source == 1 %STOXX full

    [IBassets, IBliabilities, ~, ~, ~, equity] = import_stoxx(num2str(varargin{1}));
    n=length(IBassets);
    IBvolume = sum(IBassets);
    y_assets = IBassets/IBvolume;
    y_liabilities = IBliabilities./IBvolume;

    cd(num2str(varargin{1}))
    if varargin{1}==2015
        matfile(1) = dir(fullfile(['E_node_T20.mat']));
        matfile(2) = dir(fullfile(['u_T20.mat']));
        for i = 1:length(matfile) 
            eval(sprintf('%s = cluss_collate_local(matfile(i).name);', matfile(i).name(1:end-4))); % eval evaluates function in formatspec, and gives flexible variable names  
        end
    else
        load('E_node_T20.mat');
        load('u_T20.mat');
    end
    cd ..
        
    adj = ones(n);
    for ii = 1:n; adj(ii,ii) = 0; end
    [adj, ~] = fitness_weights(adj, y_assets, y_liabilities, IBvolume); %assign weights
   
end

R1 = NaN(length(adju),1);
R1end = NaN(length(adju),1);
max_eig = NaN(length(adju),1);
h = zeros(n,T+1,length(adju));
h1 = zeros(n,T+1,length(adju));

%% create leverage matrix

adjT=adj';
L0 = NaN(n);
for i = 1:n
    L0(i,:) = adjT(i,:)./equity(i);
end

%% change equity

if pol == 0 %uniform: fixed allocation
    
    x = sum(equity)*0.05;
    equity1 = equity+x/n;

elseif pol == 1 %relative increase
    
    equity1=equity*1.05;

elseif pol == 2 % E_i based, also need line 106
    
    x = sum(equity)*e1; % sum to be allocated
    
end


%% main
% tic
for l = 1:length(adju)
%     if mod(l,10)==0
%         fprintf('at lam %d out of %d\n',l,length(adju))
%     end

    ind1=int64(adju(l)*10); %assuming adju=0.1:0.1:3;
    if source < 1
        
        u1 = u_T20(:,:,ind1);
        u = zeros(size(u1)); %make cumulative u
        for t = 1:T;
            u(:,t) = sum(u1(:,1:t),2);
        end
        E_node = E_node_T20(:,ind1);
        
    elseif source == 1
        
        if varargin{1} == 2015
            u1 = u_T20(:,:,length(xT)*(ind1-1)+ind0);
            u = zeros(size(u1)); %make cumulative u
            for t = 1:T;
                u(:,t) = sum(u1(:,1:t),2);
            end
            E_node = E_node_T20(((ind1-1)*n)+1:ind1*n,ind0);
        else
            u1 = u_T20(:,:,ind1,ind0);
            E_node = E_node_T20(:,ind1,ind0);
            u = zeros(size(u1)); %make cumulative u
            for t = 1:T;
                u(:,t) = sum(u1(:,1:t),2);
            end
        end
                
    end

    if pol == 2 %E_i based policy
      equity1 = equity+(x.*(E_node/sum(E_node))); %allocate equity increase x as a function of E_node
    end

    %reassign original leverage matrix and equity
	L=L0/max(eig(abs(L0)))*adju(l); 
	
    %calc final R1 for base   
    for t = 1:T 
        h(:,t+1,l) = L*h(:,t,l)+u(:,t);
    end
    R1end(l) = sum(h(:,end,l));
    hT(l)=min(h(:,end,l));
    
    %calc new L
    for i = 1:n
        L(i,:) = adjT(i,:)./equity1(i);
    end
    b = adju(l)/max(eig(L0));
    max_eig(l)=max(eig(abs(b*L)));
        
    %calc final R and express as fraciton of original R 
    for t = 1:T
        h1(:,t+1,l) = b*L*h1(:,t,l)+u(:,t);
    end
    R1(l) = sum(h1(:,end,l))/R1end(l);
end        
% toc


%% plotting
if p == 1
    
figure,
    p=plot(adju,R1);
    xlabel('Eigenvalue \lambda_{max}')
    ylabel(['$$ \frac{R}{R_0}$$' ' Fraction of total loss'],'interpreter','latex')

end

end