function [z, fval, exitflag] = fitness_model_param(param, y_out, y_in, k_tot)
    %%fitness model with ME, after Musmeci2012 "Bootstrapping topology and
    %%systemic risk of complex network using the fitness model"
    %%find free parameter z using fminsearch deterministic optimisation
        %linear programme, assuming global minimum exists. 
        %testing each value $param$ to minimise $err$ in function
        %fitness_model. minimum of $param$ is $z$, the output of
        %fminsearch, and the minimum of $err$ is $fval$.
        
    %%inputs
    % param = random start value, if have an ide aof likely solution,
        % speeds up process
    % y = fitness values
    % k_tot = degree of each node

    %% run model
    [z, fval, exitflag] = fminsearch(@fitness_model,param); 
    
    function err = fitness_model(param)
        
        calc = 0;
        calc_k = 0;
        for i = 1:length(y_out)
            for j = 1:length(y_out)
        
                if i~=j
                
                    calc = calc + (param.*y_out(i).*y_in(j))./(1+(param.*y_out(i).*y_in(j)));
                                                
                end
            end
        end
        
        err = sum((2*calc-sum(k_tot)).^2); % can change k

    end

end