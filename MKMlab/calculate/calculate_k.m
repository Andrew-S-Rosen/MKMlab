function k = calculate_k(A,Ea,T)
%calculate rate constant from Arrhenius equation
%
%INPUTS:
%A - double: pre-exponential factor
%Ea - double: activation energy
%T - double: absolute temperature
%
%OUTPUTS:
%k - double: rate constant

k = A*exp(-Ea/(const.R*T));

end