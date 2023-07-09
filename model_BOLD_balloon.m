function [BOLD_activity] = model_BOLD_balloon(input_ts)
% model_BOLD_balloon.m
%
% Simulate BOLD balloon model 


    nT = size(input_ts,1);
    N = size(input_ts,2);

    param = loadParameters_balloon_func;

    % Time parameters
    param.tstep    = 0.01;       % time step
    param.tmax     = nT*0.01;        % maximum time
    param.tspan    = [0, param.tmax];             % time period limits
    param.T        = 0:param.tstep:param.tmax;    % time vector

    % initialize activity vectors
    sol.z = zeros(N, length(param.T));
    sol.f = zeros(N, length(param.T));
    sol.v = zeros(N, length(param.T));
    sol.q = zeros(N, length(param.T));
    sol.BOLD = zeros(N, length(param.T));
    
    % setting initial condition
    F0 = repmat(0.001*ones(N,1), 1, 4);
    
    F = F0;
    sol.z(:,1) = F(:,1);
    sol.f(:,1) = F(:,2);
    sol.v(:,1) = F(:,3);
    sol.q(:,1) = F(:,4);
    
    for k = 2:length(param.T)
        dF = balloon_ODE(input_ts(k-1,:).', F, param);
        F = F + dF*param.tstep;
        sol.z(:,k) = F(:,1);
        sol.f(:,k) = F(:,2);
        sol.v(:,k) = F(:,3);
        sol.q(:,k) = F(:,4);
    end
    
    sol.BOLD = 100*param.V0*(param.k1*(1 - sol.q) + param.k2*(1 - sol.q./sol.v) + param.k3*(1 - sol.v));
    
    BOLD_activity = real(sol.BOLD);

end

function dF = balloon_ODE(S, F, param)
% balloon_ODE.m
%
% Calculate the temporal activity by solving an ODE.
%
% Inputs: S          : spatiotemporal external input [V x 1]
%         F          : solutions at one time point
%         param      : model parameters (struct)
%
% Output: dF         : time derivative of variables [4 x 1]
%
% Original: James Pang, Monash University, 2022

z = F(:,1);
f = F(:,2);
v = F(:,3);
q = F(:,4);

dF(:,1) = S - param.kappa*z - param.gamma*(f - 1);
dF(:,2) = z;
dF(:,3) = (1/param.tau)*(f - v.^(1/param.alpha));
dF(:,4) = (1/param.tau)*((f/param.rho).*(1 - (1 - param.rho).^(1./f)) - q.*v.^(1/param.alpha - 1));


end

function out = balloon_Fourier(input_ts, T, param)
% balloon_Fourier.m
%
% Calculate the temporal activity of one mode via Fourier transform.
%
% Inputs: mode_coeff : coefficient of the mode [1 x T]
%         T          : time vector with zero center [1 x T]
%         param      : model parameters (struct)
%
% Output: out        : activity [1 x T]
%
% Original: James Pang, Monash University, 2022

Nt = length(T);
Nw = Nt;
wsamp = 1/mean(param.tstep)*2*pi;
jvec = 0:Nw-1;
w = (wsamp)*1/Nw*(jvec - Nw/2);

% mode_coeff_fft = ctfft(mode_coeff, param.T);	
mode_coeff_fft = coord2freq_1D(mode_coeff, w);
	
% Transfer functions from Robinson et al. 2006, Aquino et al. 2012, 2014,
% Pang et al. 2016, 2018
T_Fz = 1 ./ (-(w + 1i*0.5*param.kappa).^2 + param.w_f^2);
T_yF = param.V_0 * (param.alpha*(param.k2 + param.k3)*(1 - 1i*param.tau*w) - (param.k1 + param.k2)*(param.alpha + param.beta - 1 - 1i*param.tau*param.alpha*param.beta*w))./((1 - 1i*param.tau*w).*(1 - 1i*param.tau*param.alpha*w));
T_yz = T_yF.*T_Fz;

out_fft = T_yz.*mode_coeff_fft;

% calculate inverse Fourier transform
out = real(freq2coord_1D(out_fft, w));

end



