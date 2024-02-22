function [gamma, ht_t34_predict, dhdt_t34, dhdt_t34_fitting] = ...
   calculateGamma(t3, t4, rain_mm, ht_m, lambda, hb, h_min_1, h_min_2)
   t34 = double((t3:t4) - t3); t34 = t34(:);
   ht_t34 = ht_m(t3:t4); ht_t34 = ht_t34(:);   
   rain_t34 = rain_mm(t3:t4); rain_t34 = rain_t34(:);
   % calculate derivative (dh/dt) (diffy is dh/dt of selected time range)
   dhdt_t34_ = diff(ht_t34); dhdt_t34_ = dhdt_t34_(:);
   t34_diff_ = linspace(.5, double(t4 - t3) -.5, size(dhdt_t34_(:),1)); t34_diff_ = t34_diff_(:);
   dhdt_t34 = interp1(t34_diff_, dhdt_t34_, t34, 'linear', 'extrap');
   % define cost function for solving gamma 
   % That is, finding best gamma so that, 
   % dhdt = gamma * rain - lambda * sqrt(h - hb) 
   eq28_ver = 3;
   if eq28_ver == 1
      lambda_lost_term =  lambda * sqrt(max(0., ht_t34 - hb)); 
   elseif eq28_ver == 2
      lambda_lost_term =  lambda * sqrt(max(0., ht_t34 - hb)) .* heaviside(ht_t34 - 0.3); 
   elseif eq28_ver == 3
      lambda_lost_term =  lambda * sqrt(max(0., ht_t34 - hb)) .* piecewise_linear(h_min_1, h_min_2, ht_t34); 
   end
   err_dhdt = @(gamma) ... 
		   (dhdt_t34 - (gamma * rain_t34 - lambda_lost_term));
 
   % run nonlinear least squares
   gamma_init = 0.0;
   gamma_lb = 0.0; 
   gamma_ub = 1e6;
   [gamma,resnorm,residual,exitflag,output,lsq_lambda,jacobian] = ...
       lsqnonlin(err_dhdt, gamma_init, gamma_lb, gamma_ub);
   % for debug: calculate residual 	   
%   residual_ = dhdt_t34 - (gamma * rain_t34 - lambda_lost_term);
%   norm_residual = norm(dhdt_t34 - (gamma * rain_t34 - lambda_lost_term));
%   figure; plot(t34 + double(t3), dhdt_t34, t34 + double(t3), (gamma * rain_t34 - lambda_lost_term)); legend('Measured', 'Fitting');
   % calculate predicted htdt_t34 (htdt_t34_predict) based on lambda, hb, gamma
   dhdt_t34_fitting = gamma * rain_t34 - lambda_lost_term;
   % calculate predicted ht_t34 (ht_t34_predict) based on lambda, hb, gamma
   ht_t34_predict = zeros(size(ht_t34));
   ht_t34_predict(1) = ht_t34(1);
   for i = 2: size(ht_t34_predict(:), 1)
       % predict water level step by step
       ht_t34_predict(i) = ht_t34_predict(i - 1) + ...
           gamma * rain_t34(i - 1) - lambda_lost_term(i - 1);
       % minimum of water level is zero
       if (ht_t34_predict(i) <= max(0., hb))
           ht_t34_predict(i) = max(0., hb); 
       end
   end
end


	