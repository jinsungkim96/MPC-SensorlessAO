% Function representing Von Karman spectrum as introduced in Eq. (?)
function RHP =Von_Karman_Spectr(KK)
global alfa kap_m K0 B_al
    RHP=B_al*exp(-(KK./kap_m).^2).*(KK.^2+K0^2).^(-(2+alfa)/2);
end 

