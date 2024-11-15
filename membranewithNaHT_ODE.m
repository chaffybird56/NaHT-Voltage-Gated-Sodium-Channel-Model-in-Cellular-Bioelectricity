function dYdt = membranewithNaHT_ODE(t, Y, I_parameters)
    % Unpack the gating variables
    Vm = Y(1); % Membrane potential
    m = Y(2); % Sodium activation
    h = Y(3); % Sodium inactivation
    n = Y(4); % Potassium activation
    w = Y(5); % NaHT activation
    z = Y(6); % NaHT inactivation

    %Determine current based on input parameters
    %---------------------------------------------------------------------%
    I_1 = I_parameters(1);  % Pulse 1 strength [μA/cm^2]
    L_1 = I_parameters(2);  % Pulse 1 duration [ms]
    I_2 = I_parameters(3);  % Pulse 2 strength [μA/cm^2]
    L_2 = I_parameters(4);  % Pulse 2 duration [ms]
    gap = I_parameters(5);  % Gap width [ms]
    
    I = 0;  % Initialize current to 0
    if (t>=0 && t<=L_1)  
        I = I_1;  % Set current to pulse 1 strength if within range 
    elseif (t>=(L_1+gap) && t<=(L_1+gap+L_2)  && L_2 ~= 0)
        I = I_2;  % Set current to pulse 2 strength if within range
    end
    %---------------------------------------------------------------------%

    %Intracellular ion concentration vector in [mM]
    %---------------------------------------------------------------------%
    C_in = [184.9 0.1 5 10]'; 
    %---------------------------------------------------------------------%
    % C_in(1) = Potassium intracellular concentration, CK_in
    % C_in(2) = Calcium intracellular concentration, CCa_in
    % C_in(3) = Chlorine intracellular concentration, CCl_in
    % C_in(4) = Sodium intracellular concentration, CNa_in
    %---------------------------------------------------------------------%
    
    
    %Extracellular ion concentration vector in [mM]
    %---------------------------------------------------------------------%
    C_out = [9 0.335 83.95 307.3]';
    %---------------------------------------------------------------------%
    % C_out(1) = Potassium extracellular concentration, CK_out
    % C_out(2) = Calcium extracellular concentration, CCa_out
    % C_out(3) = Chlorine extracellular concentration, CCl_out
    % C_out(4) = Sodium extracellular concentration, CNa_out
    %---------------------------------------------------------------------%
    
    
    %Ion resting conductance vector in [mS/cm^2]
    %
    %---------------------------------------------------------------------%
    g = [0.4 0.05 0.5 0.05];
    %---------------------------------------------------------------------%
    % g(1) = Potassium resting conductance, gK
    % g(2) = Calcium resting conductance, gCa
    % g(3) = Chlorine resting conductance, gCl
    % g(4) = Sodium resting conductance, gNa

    %Calculate voltage-gated conductances [mS/cm^2]
    g_Na_v = 240;  
    g_K_v = 120;
    g_Na = g_Na_v*m^3*h;  %voltage-gated sodium conductance
    g_K = g_K_v*n^4;  %voltage-gated potassium conductance
    
    %Add voltage-gated conductances to passive conductances
    g(1) = g(1) + g_K;
    g(4) = g(4) + g_Na;

    %---------------------------------------------------------------------%
    
    %Ion valence vector [unitless]
    %---------------------------------------------------------------------%
    Z = [1 2 -1 1]; %Vector of ion valences
    %---------------------------------------------------------------------%
    % Z(1) = Potassium ion valence, ZK
    % Z(2) = Calcium ion valence, ZCa
    % Z(3) = Chlorine ion valence, ZCl
    % Z(4) = Sodium ion valence, ZNa
    %---------------------------------------------------------------------%
    
    
    %Initialize equilibrium potential vector in [mV]
    %---------------------------------------------------------------------%
    E = zeros(4,1);
    %---------------------------------------------------------------------%
    % E(1) = Potassium equilibrium potential, Ek
    % E(2) = Calcium equilibrium potential, ECa
    % E(3) = Chlorine equilibrium potential, ECl
    % E(4) = Sodium equilibrium potential, ENa
    %---------------------------------------------------------------------%
    
    
    %Other variables
    %---------------------------------------------------------------------%
    R = 8.314;       %Gas constant in [J/(K*mol)]
    F = 96487;       %Faraday's constant in [C/mol]
    T = 15 + 273.15; %Membrane temperature in [K]
    C_m = 1;          %Membrane capacitance in [uF/cm^2]
    %---------------------------------------------------------------------%

    %Calculate Nernst equilibrium potentials for all ion species
    %Calculating values here allows us to change membrane parameters
    %without needing to manually do a calculation for Ei
    for i = 1:length(E)
        E(i) = (-R*T/(Z(i)*F))*log(C_in(i)/C_out(i)); %Calc Ei
        E(i) = E(i)*1000;                             %Convert to mV
    end

    E_NaHT = E(4); % Since E_NaHT = E_Na
    vm = Vm + 60;  %calculate relative membrane potential [mV]

    %---------------------------------------------------------------------%
    % Calculating particle activiation and inactiviation  
    % transition rates. These equations were derived from Hodgkin-Huxel's
    % model.

    %Sodium activation transition rates
    alpha_m = 0.374*(vm-25.41)/(1-exp((25.41-vm)/6.06));
    beta_m = 0.795*(21-vm)/(1-exp((vm-21)/9.41));

    %Sodium inactivation transition rates    
    alpha_h = -0.110*(27.74+vm)/(1-exp((27.74+vm)/9.06));
    beta_h = 4.514/(1+exp((56-vm)/12.5));

    %Potassium activation transition rates    
    alpha_n = 0.0516*(vm-35)/(1-exp((35-vm)/10));
    beta_n = 0.129*(35-vm)/(1-exp((vm-35)/10));

    
    %---------------------------------------------------------------------%
    % Define additional parameters for NaHT
    g_NaHT_bar = I_parameters(6); % Maximal conductance for high-threshold Na+ channels in mS/cm^2

    % NaHT activation and inactivation dynamics
    alpha_w = 0.0936 * (vm - 55.41) / (1 - exp((55.41 - vm) / 6.06));
    beta_w = 0.199 * (51 - vm) / (1 - exp((vm - 51) / 9.41));
    alpha_z = -0.055 * (27.74 + vm) / (1 - exp((vm + 27.74) / 9.06));
    beta_z = 2.257 / (1 + exp((56 - vm) / 12.5));

    % Compute the high-threshold Na+ current
    I_NaHT = g_NaHT_bar * w * z * (Vm - E_NaHT);
    

    %---------------------------------------------------------------------%
   
    % Ensure g is a column vector
    g = g(:);
        
        % Calculate passive currents for all ions
    total_passive_current = sum(g .* (Vm - E));
    
    % Include I_NaHT in the total membrane current calculation
    dVmdt = (I - total_passive_current - I_NaHT) / C_m;
    
    %---------------------------------------------------------------------%
    % Calculating the derivatives obtained from the given equations 
    dmdt = alpha_m*(1-m)-beta_m*m;
    dhdt = alpha_h*(1-h)-beta_h*h;
    dndt = alpha_n*(1-n)-beta_n*n;

    % Calculate derivatives for w and z
    dwdt = alpha_w * (1 - w) - beta_w * w;
    dzdt = alpha_z * (1 - z) - beta_z * z;

    % Construct the derivative vector
    dYdt = [dVmdt; dmdt; dhdt; dndt; dwdt; dzdt];
end
