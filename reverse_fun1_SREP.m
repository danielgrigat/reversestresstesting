function [E, E_node, IPR, loss_frac, h, adj, u, flag, varargout] = reverse_fun1_SREP(adju, xT, T, so, p, varargin)
%% function of reverse stress testing framework
%inputs
% adju = desired \lambda_max, e.g. single entry adju = 0.5; or vector, e.g. adju = 0.1:0.1:2;
% xT = target state, individual double, assumed that target same for all nodes
% T = control time, individual integer
% so = source:
%                 so<1 for density variable; here varargin{1} is year as string ('2014','2015','2016')
%                 so=2 for full STOXX; here varargin{1} is year as string ('2014','2015','2016')
%                 so=6 for external adj and equity via varargin
% p = plots on p=1; off p=0
% varargin{1} = external adjacency matrix || for so <1 & 2 it is year as string
% varargin{2} = external beta to adjust \lambda_max, if empty will compute itself
    %beta vector must be same size as adju, as I am looping over adju
% varargin{3} = external equity

%Output:
%varargout: vector of \lambda_max when fed external beta via varargin{2}


%N.B. quadprog and fmincon are identical, but the former is much faster

%% check inputs
if length(varargin)>1 && ~isempty(varargin{2})  
    if length(varargin{2}) ~= length(adju)
        error('Vector of betas must be same length as vector of desired lambdas')
    end
end

%% get adjacency matrix

if so < 1

    [IBassets, IBliabilities, ~, ~, ~, equity] = import_stoxx1(varargin{1});
    n=length(IBassets);
    IBvolume = sum(IBassets);
    y_assets = IBassets/IBvolume;
    y_liabilities = IBliabilities./IBvolume;
    
    d1 = so*2; %density, needs to be set as 2x desired density
    ed = round(d1*n^2); %desired number of edges
    param = ed; % starting value for free parameter estimation
    [z, fval, exitflag_p] = fitness_model_param(param, y_assets, y_liabilities, ed); %get free parameter
    adj = fitness(n, y_assets, y_liabilities, z); %get binary
    for ii = 1:n; adj(ii,ii) = 0; end
    [adj, exitflag_w] = fitness_weights(adj, y_assets, y_liabilities, IBvolume); %assign weights
    fprintf('density I am working with now is %1.4f\n',nnz(adj)/n^2);
    fprintf('fraction of weight allocated %1.4f\n',sum(adj(:))/IBvolume);
    
elseif so == 2

    [IBassets, IBliabilities, ~, ~, ~, equity] = import_stoxx1(varargin{1});
    n=length(IBassets);
    IBvolume = sum(IBassets);
    y_assets = IBassets/IBvolume;
    y_liabilities = IBliabilities./IBvolume;    
    
    adj = ones(n);
    for ii = 1:n; adj(ii,ii) = 0; end
    [adj, exitflag_w] = fitness_weights(adj, y_assets, y_liabilities, IBvolume); %assign weights

    
elseif so == 6 && ~isempty(varargin)
    
    adj = varargin{1};
    equity = varargin{3};
    n = length(equity);
    
end


%% set up leverage matrix

adjT=adj';
L = NaN(n); %leverage matrix, %makes no difference to include active banks as in paper.
for i = 1:n
    L(i,:) = adjT(i,:)./equity(i);
end


%% control


Beq=ones(n,1).*xT; %fmincon, target state
x0=zeros(n*T,1);
options = optimoptions('quadprog','tolfun',10^-15,'TolCon',10^-15,'MaxIterations',10000000,'Display','off');

u = NaN(n,T,length(adju));
IPR = NaN(length(adju),1);
loss_frac = zeros(length(adju),1);
h=zeros(n,T+1,length(adju));
E_node = NaN(n,length(adju));

for i = 1:length(adju)
    if mod(i,10)==0
        fprintf('at lambda %d out of %d\n',i,length(adju))
    end

    if length(varargin)>1 && ~isempty(varargin{2})
        LL = L*varargin{2}(i);  
        varargout{i,1} = max(eig(abs(LL)));
    else
        b = adju(i)/max(eig(abs(L)));
        LL = L*b;        
    end
    %%%%%%%%%%%%%
%     fprintf('average edge weight: %10.10f\n', mean(LL(:)))
%     fprintf('volume: %10.10f\n', sum(LL(:)))
%     fprintf('volume: %10.10f\n', sum(LL(:)))
    %%%%%%%%%%%%%
    Aeq=[];    
    for s=1:T
        L1=0;
        for t=s:T
            L1 = L1 + LL^(T-t);
        end
        Aeq = [Aeq L1];
    end

%     [y,E(i),flag(i), output(i)] = fmincon(fun,x0,-Aeq,-Beq,[],[],[],[],[],options); %reach at least certain level Beq
%     [y,E(i),flag] = fmincon(fun,x0,[],[],Aeq,Beq,[],[],[],options); %equality: Aeq*u = Beq    
    [y,E(i),flag(i), output(i)] = quadprog(eye(n*T)*2,x0,-Aeq,-Beq,[],[],[],[],[],options); %reach at least certain level Beq
    
	u(:,:,i)=reshape(y,[n,T]);
	E_node(:,i) = sum(u(:,:,i).^2,2);

    for t = 1:T;
        u1 = sum(u(:,1:t,i),2); %cumulative u
        h(:,t+1,i) = LL*h(:,t,i)+u1;
    end    
    
    
    IPR(i)=1/sum((E_node(:,i)/E(i)).^2);
    loss_frac(i) = sum(sum(u(:,:,i),2).*equity)/sum(equity);

end

%% plotting

if p == 1

figure, plot(adju,E)
xlabel('Largest eigenvalue \lambda_{max}')
ylabel('Energy E')

figure, plot(adju,IPR)
xlabel('Largest eigenvalue \lambda_{max}')
ylabel('IPR')

figure, plot(adju,loss_frac)
xlabel('Largest eigenvalue \lambda_{max}')
ylabel('Loss frac')

end

end
