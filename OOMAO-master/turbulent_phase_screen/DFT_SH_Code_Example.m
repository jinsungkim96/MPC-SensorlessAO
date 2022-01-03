%Example of Matlab code generating N_MC phase samples 
%at N_FT by N_FT square grid using DFT and sub-harmonics correction
%Von Karman spectrum with inner scale l0, and outer scale L0 is used.
%DFT with N_FFT*N_FFT terms and SH correction of N_SH order.
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
%DFT wavenumbers
k=(-N_FT/2:N_FT/2-1)*dk;
K=sqrt(repmat(k,N_FT,1).^2+repmat(k',1,N_FT).^2);
%RMS values of spectral amplitudes
Sig_DFT=dk*sqrt(Von_Karman_Spectr(K));Sig_DFT(K==0)=0;
%-------------- Sub Harmonics Argument arrays -----------------------
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
%3 by 3 by N_SH arrays of Cartesian components of SH wavevectors   
    k_SH_x=dk_SH_3.*SH_Gr_x_3;
    k_SH_y=dk_SH_3.*SH_Gr_y_3;
%9*N_SH vector of Cartesian components of SH wavevectors      
    k_SH_x1=reshape(k_SH_x,9*N_SH,1);
    k_SH_y1=reshape(k_SH_y,9*N_SH,1)';
%9*N_SH vector of SH wavenumbers      
    k_SH=sqrt(k_SH_x1.*k_SH_x1+(k_SH_y1.*k_SH_y1)');
%9*N_SH vector of RMS values of SH amplitudes      
    Sig_SH1=dk_SH_1.*sqrt(Von_Karman_Spectr(k_SH));
%Remove SH overlaps and DC contribution
    Sig_SH1(5:9:9*N_SH)=0;     %No zero contribution
%Auxilliary  9*N_SH by N_x and N_x by 9*N_SH matrices used for SH phase 
%calculations inside the Monte-Carlo loop
    Exp_x1=exp(1i*k_SH_x1*x); Exp_y1=exp(1i*y*k_SH_y1);
end  
%------------------------ MC Loop ------------------------------
tic
for i_MC=1:N_MC
% - - - - - - - DFT part of phase- - - - - - -
%Random complex amplitude of DFT spectral components
	Ampl_DFT=Sig_DFT.*(randn(N_FT)+1i*randn(N_FT));
%Complex and real DFT phases as N_FT by N_FT matrices  
    psi_DFT=fftshift(fft2(fftshift(Ampl_DFT)));
    fi_DFT1=real(psi_DFT);fi_DFT2=imag(psi_DFT);
%- - - - - - - Subharmonics part of phase - - - - - - -
    if N_SH>0 
%9*N_SH vector of random complex SH amplitudes         
        Amp_SH=Sig_SH1.*(randn(9*N_SH,1)+1i*randn(9*N_SH,1));
%N_FT by N_FT matrix of complex SH phase       
        psi_SH=Exp_y1*diag(Amp_SH)*Exp_x1;
    else
        psi_SH=zeros(N_FT); 
    end
%Two N_FT by N_FT matrices of real SH phases 
    fi_SH1=real(psi_SH); fi_SH2=imag(psi_SH);
%Two N_FT by N_FT matrices of real DFT-SH phases 
    fi1=fi_DFT1+fi_SH1;fi2=fi_DFT2+fi_SH2;
end %MC
toc

