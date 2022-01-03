%Example of Matlab code generating N_MC phase samples 
%at N_x by N_x square grid using SS technique with N_SS spectral components
%Von Karman spectrum with inner scale l0, and outer scale L0 is used.
%Log-uniform partition of wave numbers of the [K_Min, K_Max] range.
%Uniform distributions of wave vectors in individual "rings," as discussed
%in Section 2.2
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
   
%-------------- SS phase parameters -----------------------  
N_SS=500; %Number of SS spectral components used
%Log-normal partition of spectral range [K_Min,K_Max]
if K_Min==0
%For K_Min=0 use K0/15 as lowest non-zero wavenumber     
    k_Bdr=[0 exp(linspace(log(K0/15),log(K_Max),N_SS))]';
else
    k_Bdr=[exp(linspace(log(K_Min),log(K_Max),N_SS+1))]';
end
%Integrand for SS spectral weights 
F_Sig=@(t)t.*Von_Karman_Spectr(t); 
%RMS values of spectral amplitudes
Sig_SS=zeros(N_SS,1);
for i=1:N_SS
    Sig_SS(i)=sqrt(2*pi*integral(F_Sig,k_Bdr(i),k_Bdr(i+1),'RelTol',1.0e-3));
end
tic
%------------------------ MC Loop ------------------------------
for i_MC=1:N_MC
%Uniform wavevectors in the annulae
    k_SS=sqrt(k_Bdr(1:N_SS).^2+(k_Bdr(2:N_SS+1).^2-k_Bdr(1:N_SS).^2).*rand(N_SS,1));
    tet_SS=rand(N_SS,1)*2*pi;
%Cartesian coordinates of wavevectors    
    kx_SS=k_SS.*cos(tet_SS); ky_SS=(k_SS.*sin(tet_SS))';
%Random complex amplitudes of spectral components    
    Amp_SS=randn(N_SS,2)*[1;1i].*Sig_SS;  
%Sample of complex SS phase
    psi_SS=exp(1i*x'*ky_SS)*diag(Amp_SS)*exp(1i*kx_SS*x);
%Two samples of real SS phase    
    fi_SS1=real(psi_SS);  fi_SS2=imag(psi_SS);
end
toc


  