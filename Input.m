global L E_s gama_p w_l w_sdl fy beta1 gama_conc
% Beam span(m)
L = 40;
% Modulus of elasticity of steel (Mpa)
E_s = 200000;
% Factor for the type of prestressing tendon
gama_p = 0.4;
% Unfactored live load per unit length (kN/m^2)
w_l = 1.75;
% Dead load due to superimposed loading plus sustained live load (kN/m^2)
w_sdl = 0.5;
% Specified yield strength of reinforcement (MPa)
fy = 400;
% Betonun ozgur agirligi (kN/m^3)
gama_conc = 25;
% Factor relating depth of equivalent rectangular compressive stress block to neutral axis depth
beta1 = 0.85;