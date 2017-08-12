function [ adj, exitflag, er_even, er_odd ] = fitness_weights(adj, y_assets, y_liabilities, totalInterbankVolume)
%% RAS - assign weights fitness model, after Battiston 2015 / Gietl and Reffel 2013, Pukelsheim 2013, Bacharach 1965
%convergence criterion for odd step based on even step, vice versa.

%%INPUTS
% adj = binary adjacency matrix found via fitness.m
% y_assets = fitness vector for all assets
% y_liabilities = fitness vector for all liabilities
% totalInterbankVolume = totalInterbankVolume of real network

%%OUTPUTS
% adj = estimate of weighted directed adjacency matrix
% exitflag = convergence flag, 1 for convergence, 0 for non convergence
% er_even = vector of error scores at each iteration for each node for columns
% er_odd = vector of error scores at each iteration for each node for rows


%% run model

iteration = 0;      
max_iteration = 1000;

er_odd = zeros(max_iteration,length(adj)); % error for odd step
er_even = zeros(max_iteration,length(adj)); % error for even step
conv_crit_e(1:nnz(adj),1)=1;
conv_crit_o(1:nnz(adj),1)=1;
conv_upper_o = log10(max(y_assets)); %initial criterion for convergence
cond_break = ones(2,1);

while sum(cond_break) ~= 0 && iteration < max_iteration

    iteration = iteration + 1;

    if mod(iteration,10) == 0
    
        fprintf(['Iteration: %d ...\n'],iteration); 
        
    end
        
        
	%% step 1 (odd step / y_liabilities / rows)
    
    for i = 1:length(adj)
        if sum(adj(i,:)) > 0
            adj(i,:) = (adj(i,:)/sum(adj(i,:)))*y_liabilities(i);
            er_odd(iteration,i) = sum(adj(i,:))-y_liabilities(i);
        end                       
                  
    end    
   
    %convergence test
    cond_break(1) = sum(er_odd(iteration,:)) > conv_upper_o;
    %fprintf(['conv_upper_o: %d\n'],conv_upper_o); 
    %fprintf(['abs error: %d\n'],abs(sum(er_odd(iteration,:)))); 
    %fprintf(['cond break(1): %d\n'],cond_break(1));        
    
    %convergence criterion for even step (y_assets)
    conv_upper_e = abs(log10(max(adj(adj>0)./conv_crit_e(conv_crit_e>0))));
    conv_crit_e = adj;
                
    
    %% step 2 (even step / y_assets / columns)
    
    for i = 1:length(adj)
        
        if sum(adj(:,i)) > 0
            
            adj(:,i) = (adj(:,i)/sum(adj(:,i)))*y_assets(i);
            er_even(iteration,i) = sum(adj(:,i))-y_assets(i); 
            
        end             
                    
    end
    
    %convergence test
	cond_break(2) = sum(er_even(iteration,:)) > conv_upper_e;
	%fprintf(['conv_upper_e: %d\n'],conv_upper_e); 
    %fprintf(['abs error: %d\n'],abs(sum(er_even(iteration,:)))); 
    %fprintf(['cond break(2): %d\n'],cond_break(2));  
    
    %convergence criterion for odd step (y_liabilities)
    conv_upper_o = abs(log10(max(adj(adj>0)./conv_crit_o(conv_crit_o>0))));
    conv_crit_o = adj;
           
        
end

%assign volume to relative matrix
adj = adj*totalInterbankVolume;
    
% exitflag
if iteration < max_iteration
    exitflag = 1; % converged
else
    exitflag = 0; % max iterations reached

end   

end
