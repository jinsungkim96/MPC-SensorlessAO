%Example of Matlab code generating N_MC phase samples 
%at N_FT by N_FT square grid using PWD and sub-harmonics correction
%Von Karman spectrum with inner scale l0, and outer scale L0 is used.
%PWD with N_FFT*N_FFT terms and SH correction of N_SH order.
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
%-------------- DFT Arguments -----------------------
N_FT=1024; %Number of spectral components used
N_SH=4; %Number of subharmonic layers
x_Max=1.0;
dx=x_Max/N_FT; 
dk=2*pi/N_FT/dx; 
K_Max=dk*N_FT/2;K_Min=dk;
x=(-N_FT/2:N_FT/2-1)*dx;
y=x';
k=(-N_FT/2:N_FT/2-1)*dk;
%-------------- Auxilliary Sub Harmonics Argument arrays -----------------------
if N_SH>0
%Create 9*N_SH vector of SH spectral cells sizes       
    dk_SH=dk./3.^[1:N_SH];
    dk_SH_3=reshape(repmat(dk_SH,9,1),3,3,N_SH);
    dk_SH_1=reshape(dk_SH_3,9*N_SH,1);
%Auxulliary 3-D grids
    SH_Gr_x=repmat([-1 0 1],3,1);
    SH_Gr_y=SH_Gr_x';
    SH_Gr_x_3=repmat(SH_Gr_x,1,1,N_SH);
    SH_Gr_y_3=repmat(SH_Gr_y,1,1,N_SH);
end
%------------------------ MC Loop ------------------------------
tic
for i_MC=1:N_MC
% - - - - - - - PWD part - - - - - - -
% Randomly displaced components of DFT wavevectors kx and ky
    kx_r=dk*(rand-0.5);
    ky_r=dk*(rand-0.5);
    kx=k+kx_r;ky=k'+ky_r;
% N_FT by N_FT correction factor of PWD technique    
    E_PWD=exp(-1i*ky_r.*y)*exp(-1i*kx_r.*x);
% Random RMS amplitudes of DFT spectral components    
    K2=repmat(kx,N_FT,1).^2+repmat(ky,1,N_FT).^2;
    Sig_PWD=dk*sqrt(Von_Karman_Spectr(sqrt(K2)));
%Avoid double account for zero area, but keep zero term for pure PWD
    if N_SH>0 
        Sig_PWD(N_FT/2+1,N_FT/2+1)=0;
    end
%Random complex amplitude of PWD spectral components    
    Ampl_PWD=Sig_PWD.*(randn(N_FT)+1i*randn(N_FT));
%Complex and real PWD phases as N_FT by N_FT matrices      
    psi_PWD=fftshift(fft2(fftshift(Ampl_PWD)));
    psi_PWD=E_PWD.*psi_PWD;
    fi_PWD1=real(psi_PWD);fi_PWD2=imag(psi_PWD);
     %******************* Subharmonics ************************
    if N_SH>0 
%3 by 3 by N_SH arrays of Cartesian components of random PWD wavevectors          
        k_SH_x=dk_SH_3.*(SH_Gr_x_3+(rand(3,3,N_SH)-0.5));
        k_SH_y=dk_SH_3.*(SH_Gr_y_3+(rand(3,3,N_SH)-0.5));
%9*N_SH vector of Cartesian components of random PWD wavevectors    
        k_SH_x1=reshape(k_SH_x,9*N_SH,1);
        k_SH_y1=reshape(k_SH_y,9*N_SH,1)';
%9*N_SH vector of RMS values of SH spectral amplitudes    
        k_SH_2=k_SH_x1.*k_SH_x1+(k_SH_y1.*k_SH_y1)';
        Sig_SH=dk_SH_1.*sqrt(Von_Karman_Spectr(sqrt(k_SH_2)));
%Avoid overlap of subharmonics, but keep zero contribution for highest subharmonic
    	Sig_SH(5:9:9*(N_SH-1))=0; 
%9*N_SH vector of random complex SH amplitudes         
        Amp_SH=Sig_SH.*(randn(9*N_SH,1)+1i*randn(9*N_SH,1));
        psi_SH=exp(1i*y*k_SH_y1)*diag(Amp_SH)*exp(1i*k_SH_x1*x);
    else   
        psi_SH=zeros(size(psi_PWD));
    end
%Two N_FT by N_FT matrices of real SH phases     
    fi_SH1=real(psi_SH);fi_SH2=imag(psi_SH); 
%Two N_FT by N_FT matrices of real PWD-SH phases 
    fi1=fi_PWD1+fi_SH1;fi2=fi_PWD2+fi_SH2;
end %MC
toc