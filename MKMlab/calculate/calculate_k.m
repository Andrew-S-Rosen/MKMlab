function k = calculate_k(A,Ea,T)
%calculate rate constant from Arrhenius equation

k = A.*exp(-Ea./(const.R*T));

end