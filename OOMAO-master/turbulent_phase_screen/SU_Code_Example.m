%Example of Matlab code generating N_MC phase samples 
%at N_x by N_x square grid using SU technique with N_SU spectral components
%Von Karman spectrum with inner scale l0, and outer scale L0 is used.
%Log-uniform partition of wave numbers of the [K_Min, K_Max] range.

%-------------- Global parameters -------------------
global alfa kap_m K0 B_al
%Global parameters used by "Von_Karman_Spectr" function
N_MC=1000;
%------------- Spectral parameters --------
alfa=5/3;    %exponent for Structure function. 1 < alfa < 2
L0=10;       %Outer scale, m
l0=0.001;    %Inner scale, m
kap_m=2*pi/l0; %High freq spectral cutoff, 1/m
K0=2*pi/L0;    %Low freq spectral cutoff, 1/m
r_C=1.0;       %coherence radius, m
% Forefactor for spectral density
B_al=r_C^(-alfa)*alfa*gamma(1+alfa/2)*2^(alfa-2)/(pi*gamma(1-alfa/2));
K_Min=K0/10;   %Low-frequency limit for SS algorithm. Can be zero.
K_Max=kap_m*2; %Hi-frequency limit for SSalgorithm.
%-------------- Spatial arguments -----------------------
N_x=1024;    %Number of spatial points in one direction
x_Max=1.0;   %Linear screen size, m
x=(1:N_x)*x_Max/N_x; %Linear array of point coordinates
%-------------- SU phase parameters -----------------------  
N_SU=500; %Number of SU spectral components used
%Log-normal partition of spectral range [K_Min,K_Max]
if K_Min==0
%For K_Min=0 use K0/15 as lowest non-zero wavenumber     
    k_Bdr=[0 exp(linspace(log(K0/15),log(K_Max),N_SU))]';
else
    k_Bdr=[exp(linspace(log(K_Min),log(K_Max),N_SU+1))]';
end
tic
%------------------------ MC Loop ------------------------------
for i_MC=1:N_MC
%Uniform wavevectors in the annulae
    k_SU=sqrt(k_Bdr(1:N_SU).^2+(k_Bdr(2:N_SU+1).^2-k_Bdr(1:N_SU).^2).*rand(N_SU,1));
    tet_SU=rand(N_SU,1)*2*pi;
%Cartesian coordinates of wavevectors    
    kx_SU=k_SU.*cos(tet_SU); ky_SU=(k_SU.*sin(tet_SU))';
%Random complex amplitudes of spectral components
    Sig_SU=sqrt(Von_Karman_Spectr(k_SU)*pi.*(k_Bdr(2:N_SU+1).^2-k_Bdr(1:N_SU).^2));
    Amp_SU=randn(N_SU,2)*[1;1i].*Sig_SU; 
%Sample of complex SS phase
    psi_SU=exp(1i*x'*ky_SU)*diag(Amp_SU)*exp(1i*kx_SU*x);
%Two samples of real SS phase    
    fi_SU_1=real(psi_SU);  fi_SU_2=imag(psi_SU);
end
toc
