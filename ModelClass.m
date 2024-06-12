classdef ModelClass
   properties
        z     = 1     ; % Instant utility of being unemployed
        delta = 0.05  ; % Probability of job destruction
        beta  = 0.5   ; % Bargain power of workers
        mu    = 5     ; % Mean of productivity
        sigma = 2.5   ; % Standard deviation of productivity
        c     = 1     ; % Cost of openning a vacancy
        rho   = 0.067 ; % Discount factor
        gamma = 0.5   ; % Matching elasticity
        lambda = 0.5  ; % Adjustment parameter
        x0=1          ; % Initial guess for fsolve
        theta_0=1     ; % Initial guess for market tightness
        dif_theta=1   ; % Initial difference for theta
        tol =  1e-6   ; % Tolerance parameter
   end
   methods
        function val = unem_solver(rhoU,z,theta,gamma,beta,mu,sigma,rho,delta)
            
            aux = @(x,beta,rhoU,rho,delta,mu,sigma) beta.*(x-rhoU)/(rho+delta)...
                                                    .*normpdf(x,mu,sigma);

            M = 1e+3;

            integral_worker = @(rhoU) integral( @(x) aux(x,obj.beta,rhoU,obj.rho,obj.delta,obj.mu,...
                                              obj.sigma),rhoU,M); % Attention with this!

            val = rhoU - obj.z - theta^(1-obj.gamma) .*  integral_worker(rhoU);
        end

        function val = theta_solver(rhoU,c,theta,gamma,beta,mu,sigma,rho,delta)

        aux = @(x,beta,rhoU,rho,delta,mu,sigma) beta.*(x-rhoU)/(rho+delta)...
                                                .*normpdf(x,mu,sigma);

        M=1e+3;

        integral_costs = @(rhoU) integral( @(x) aux(x,obj.beta,rhoU,obj.rho,obj.delta,obj.mu,...
                                          obj.sigma),rhoU,M);  % Attention with this!

        val = obj.c - theta.^(1-obj.gamma)/theta .* integral_costs(rhoU);
        end
       
       function val = run_model(obj)

        it = 0;

        options = optimoptions('fsolve', 'Display', 'off');
           
        theta_old  = obj.theta_0; %Attention with this!
        diff_theta = obj.dif_theta; %Attention with this!

        while diff_theta>obj.tol
            fun = @(rhoU) unem_solver(rhoU,obj.z,theta_old,obj.gamma,obj.beta,obj.mu,obj.sigma,obj.rho,obj.delta);
            [rhoU1, ~] = fsolve(fun,obj.x0, options);

            fun = @(theta) theta_solver(rhoU1,obj.c,theta,obj.gamma,obj.beta,obj.mu,obj.sigma,obj.rho,obj.delta);
            [theta1, ~] = fsolve(fun,theta_old, options);

            diff_theta=abs(theta_old-theta1);

            theta_old=theta_old*obj.lambda + theta1*(1-obj.lambda);
            it = it + 1;
        end

        val = [theta_old , rhoU1];
           
       end
   end
end
