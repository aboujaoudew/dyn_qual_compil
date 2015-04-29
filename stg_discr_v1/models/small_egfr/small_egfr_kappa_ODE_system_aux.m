function dydt=small_egfr_kappa_ODE_system_aux(t,y)


global perturbation_trigger;
global k;

k_temp=k;
obs_temp=small_egfr_kappa_ODE_system_obs(y);



dydt=zeros(small_egfr_kappa_ODE_system_size(),1);

