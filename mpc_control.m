%% 210930

% mpc_control_v7_v2 : Validation set 사용 + Constrained QP (via CVX)
% mpc_control_v7_v3 : Validation set 사용 + Closed-form solution (via lsqminnorm) with constraints
% mpc_control_v7_v4 : Validation set 사용 + Constrained QP (via fmincon)
% mpc_control_v7_v5 : Validation set 사용 + Constrained QP (via CasADi)

% mpc_control_v8_v3 : VAR 모델 차수 2 + Closed-form solution (via pinv) with constraints
% mpc_control_v9_v3 : DM actuator - active set 사용 (전체 144개 중 88개만 사용) -> 결정할 제어 입력의 수 감소
% dp_control_v1 : 제어 방법으로 DP를 사용

% mpc_control_v10_v2 : n = 28, m = 88 로 설정 + MPC
% mpc_control_v10_v3 : n = 28, m = 88 로 설정 + MPC

% mpc_control_v10_v3_v1 : n = 28, m = 88 로 설정 + fastMPC
% mpc_control_v10_v3_v2 : n = 28, m = 88 로 설정 + CVX (SDPT3) - 상태변수를 제약조건의 형태로 쌓지 않음
% mpc_control_v10_v3_v4 : n = 28, m = 88 로 설정 + fmincon (interior-point) - 상태변수를 제약조건의 형태로 쌓지 않음
% mpc_control_v10_v3_v5 : n = 28, m = 88 로 설정 + YALMIP (OSQP) - 상태변수를 제약조건의 형태로 쌓지 않음

% input_data_v2 : B_act matrix가 (0,1) [V] 범위로 설정된 형태
% input_data_v3 : B_act matrix가 DM spec을 고려하여 (-16,16) [rad] 범위로 설정된 형태
% input_data_v4 : B_act matrix가 DM spec을 고려하여 (-1.5,1.5) [rad] 범위로 설정된 형태
% input_data_1layer_v1 : A가 single layer phase screen (f_r0 = 1.0)을 기반으로 생성된 모델이며, B_act matrix가 DM spec을 고려하여 (-1.5,1.5) [rad] 범위로 설정된 형태
% input_data_3layer_v1 : A가 three layer phase screen (h = [5000, 5000, 5000])을 기반으로 생성된 모델이며, B_act matrix가 DM spec을 고려하여 (-1.5,1.5) [rad] 범위로 설정된 형태

% PSF에 대한 unit regularization이 진행된 데이터를 사용 (AU = 1e+12를 사용하여 [W/m^2] -> [W/um^2]로 변환)

%% Linear MPC Simulation via CVX

clear;clc;close all;

% load("./input_data/new/input_data_3layer_1000_v3.mat")
% load("./input_data/new/model_A_28.mat")

load("./input_data/new/input_data_3layer_mainly_frozen_v2.mat")
load("./data_result/mpc/mainly_frozen_flow/new_multi_n28/model_A_28_v2.mat")

load Zs_28_512.mat; % Zernike polycnomials for Estimator
load Zs_91_512.mat; % Zernike polycnomials

% load("./data_result/mpc/mainly_frozen_flow/new_multi_n28/model_A_28.mat")
load("./input_data/new/model_B_28_88.mat")
ad_acc(:,29:end) = [];

load("./data_result/SNR/SNR_10.mat")
sample = 1;

mag_conv_5 = 1;
mag_conv_7p5 = 1.401982897784687;
mag_conv_10 = 1.781797436291855;
mag_conv_12p5 = 2.145935547352314;
mag_conv_15 = 2.498049532979032;
mag_conv_20 = 3.174802103932926;
mag_conv_22p5 = 3.502222723017048;
mag_conv_25 = 3.823622456658651;

phase = phase .* mag_conv_10;
ad_acc = ad_acc .* mag_conv_10;

A = A';
A1 = A1'; % input_data에 저장된 A가 x*A 형태로 쓰였기 때문에 변환 필요
A2 = A2';
% A = eye(90);

% piston element is removed
B(1,:) = []; A_s(:,1) = []; A_s_est(:,1) = []; ad_acc(:,1) = [];
% ad_acc = ad_acc/2;

%수차 크기가 1이 되도록 normalization (자기 자신의 norm으로 나눔)
% for j=1:size(ad_acc,1)
%     ad_acc(j,:) = ad_acc(j,:) / norm(ad_acc(j,:));
% end

%% Setup the simulation parameters
nx = size(A1,1); % number of state
nu = size(B,2); % number of input
p = size(A_s,1); % number of measurements
nx_est = size(A_s_est,2);

n_mode = 6; % number of zernike modes
N = 2; % prediction horizon for MPC
T_final = 500; % 500
% T_final = size(phase,3) - num_train - num_test;
T_settle = 1.0e-3;
T_s = 0.1e-3; % sampling time of MPC Simulation
T_s_tur = 5e-3;

%% J = U'*H*U + r'*U + c
Q = (1.5e+4)*eye(nx,nx); % weighting matrix for error state
P = (1e+0)*Q; % weighting matrix for terminal error
R = (1e+0)*eye(nu,nu);% weighting matrix for cost function
% W_coeff = 0.01; % process noise -> 굳이 추가될 필요 X (이미 model ID 등으로부터 fitting error가 포함됨)
disturbance_coeff = 0.0; % 5e-1;

% [K, P, e] = dlqr(A,B,Q,R); % P is discrete-time ARE (DARE) Solution

% SNR = 10; % the magnitude of desired SNR [db]

coeff_a = 0.047275; coeff_b = 2.709264; coeff_c = 0; % Unit change parameters (Voltage to Deflection)

u_min = -28*ones(nu,1); % input box constraint [rad]
u_max = 28*ones(nu,1);  % 28 [rad] -> 200 [V]

% du_min = u_min(1,1)*(T_s/T_settle)*ones(nu,1); % input constraint
% du_max = u_max(1,1)*(T_s/T_settle)*ones(nu,1);

du_min = -0.0662*ones(nu,1); % ramp-rate constraint [rad]
du_max = 0.0662*ones(nu,1);

% du_min = -0.1739*ones(nu,1); % ramp-rate constraint [rad]
% du_max = 0.1739*ones(nu,1);

U_min = repmat(u_min,N,1);
U_max = repmat(u_max,N,1);

dU_min = repmat(du_min,N,1);
dU_max = repmat(du_max,N,1);

%% Make pixel array for image plane and pupil plane

% y_ref generation by FFT2
len = 512; %2048*2; % number of pixels in the array
mag = 1;
res = len*mag; % resolution of FFT
cen = len/2 + 1;
dx = 6.5e-6; %0.1e-6;   % pixel spacing (m)
df = 1/(len*dx); % spacing in the spatial frequency domain (cycles/m)
wavelength = 532.0e-9;   % meter
unit_change = wavelength/(2*pi)*1e+9; % [rad] to [nm]
    
xaxis = ((-len/2):(len/2-1))*dx;
yaxis = -xaxis;
xaxis_res = ((-res/2):(res/2-1))*dx/mag;
yaxis_res = -xaxis_res;

x = (-(len-1):2:(len-1))/(len-1);
[X,Y] = meshgrid(x);
[theta,r] = cart2pol(X,Y);    % convert to polar coordinate
is_in = r <= max(abs(x));
r = r(is_in);
theta = theta(is_in);

range_min = find(abs(xaxis_res-(-1.0e-4)) < 3e-6, 1, 'first');
range_max = find(abs(xaxis_res-(+1.0e-4)) < 3e-6, 1, 'last');
diff = (range_max - range_min + 1);
AU = 1e+12;

fxaxis = ((-len/2):(len/2-1))*df;
fyaxis = -fxaxis;
[FX,FY] = meshgrid(fxaxis,fyaxis);  %2-D arrays hold fx location and fy location of all points
freq_rad = sqrt(FX.^2 + FY.^2);
maxfreq = (len/2-1)*df;

pupil_radius = 1.0*maxfreq; %NA / wavelength % 1.0*maxfreq; % radius of the pupil, inverse microns % (Hanser, 2004) doi.org/10.1111/j.0022-2720.2004.01393.x
pupil_area = pi*pupil_radius^2; % pupil area
pupil = double(freq_rad <= pupil_radius); % pin-hole

idx2 = 5;

zd_dist = 3;
zd_list = (-zd_dist:zd_dist:zd_dist); %(-2:0.04:2);  %(-0.5:0.01:0.5); %distance of zd_list


%% Design matrix generation for Linear MPC

[H, M1, M2, Q_tilda, R_tilda, B_conv, B_conv_pre1, B_conv_pre2, E] = MPC_DesignMatrices_v4(A1,A2,B,nx,nu,N,Q,P,R);

closed_form_matrix = -0.5*(pinv(H'*H))*(H');

%% MPC Simulation for total simulation time

U_acc = []; U_v_acc = []; dU_acc = []; X_acc = []; X_acc_err = [];
X_err = zeros(T_final,N); X_est_err = []; X_err_low = []; X_err_high = [];
ad_est_acc = []; ad_cor_acc = []; J = [];
T_mldivide = zeros(T_final,1); T_pinv = zeros(T_final,1); T_lsqmin = zeros(T_final,1); T_geninv = zeros(T_final,1); T_cvx = zeros(T_final,1); T_fmincon = zeros(T_final,1); T_sim = zeros(T_final,1);
solver_objective = zeros(T_final,1); solver_feasibility = cell(T_final,1);


phase_cor = zeros(len,len,T_final); % corrected aberration by DM
phase_res = zeros(len,len,T_final); % residual aberration (phase_res = phase + phase_cor)
RMS_phase_res = zeros(T_final,1);

u_prev = zeros(nu,1); Y_M_acc = [];

ad_acc_valid = ad_acc(1+num_train+num_test:end,:);
phase_valid = phase(:,:,1+num_train+num_test:end);
% ad_acc_valid = ad_acc(401:800,:);
% phase_valid = phase(:,:,401:800);

for iSimStep = 1:T_final
    tic
    if iSimStep == 1
        phase_res(:,:,iSimStep) = phase_valid(:,:,iSimStep);
    else
        % ramp constraint for first control input, u[0|k]
        dU_min(1:nu,1) = du_min + u_prev;
        dU_max(1:nu,1) = du_max + u_prev;
        
        phase_res(:,:,iSimStep) = phase_valid(:,:,iSimStep) + phase_cor(:,:,iSimStep-1);
    end
    
    %% Estimator for Residual Aberration
    % Estimator (Least-square via 1st approximation model)
    Y_M = []; Y_M_noise = [];
    
    z = phase_res(:,:,iSimStep); % phase screen, unit : [m] (magnitude -1e-7 ~ -1e-7)
    % z = phase_res(:,:,iSimStep)/2; % phase screen, unit : [m] (magnitude -1e-7 ~ -1e-7)
    scrn = z;
    
    for k=1:size(zd_list,2) %zd_list = -0.5:0.01:0.5 %-0.25:0.01:0.25 %-0.5:0.005:0.5 %-1:0.05:1
        zd_diversity = zd_list(k);
        
        kW = zd_diversity.*squeeze(Zs(idx2,:,:)); % [rad] 비교를 동일하게 해주기 위해 phase diversity를 사용
        P_defocus = pupil.*exp(1i*(scrn+kW));
        
        % transfer function
        I_defocus = fftshift(fft2(ifftshift(P_defocus),res,res))*dx^2;
        im = abs(I_defocus).^2;
        v_im(:,:,k) = im(range_min:range_max,range_min:range_max)*AU;
        Y_M = [Y_M; reshape(v_im(:,:,k),(range_max-range_min+1)^2,1)];
        
        % Y_M_noise = [Y_M_noise; rand(diff^2,1)*mean(Y_M((k-1)*diff^2+1:k*diff^2,1))/SNR];
    end
    Y_M_noise = squeeze(SNR_data(sample,iSimStep,:));
    
%     Y_M = Y_M + Y_M_noise;
    Y_M_acc = [Y_M_acc; Y_M'];
    
    % ad_est = (A_s_est'*A_s_est)\(A_s_est)'*(Y_M-b_s_est); % ad_est = [0; (A_s'*A_s)\(A_s)'*(Y_M-b_s)];
    ad_est = lsqminnorm((A_s_est'*A_s_est),((A_s_est)'*(Y_M-b_s_est)));
%     if iSimStep==1 || iSimStep==2 || iSimStep==3
%         ad_est = ad_acc_valid(iSimStep,:)';
%     else
%         ad_est = lsqminnorm((A_s_est'*A_s_est),((A_s_est)'*(Y_M-b_s_est)));
%     end
    
    
    X_est_err = [X_est_err; norm(ad_est,2)];
    
    
    % residual aberration에 대한 zernike coefficient
    ad_est = [ad_est; zeros(nx-nx_est,1)]; % random value : 1e-1*2*(rand(15,1)-0.5)
    ad_est_acc = [ad_est_acc; ad_est'];
    
    % initial condition for MPC Controller
    x0 = ad_est;
    if iSimStep == 1
        x0_pre = zeros(nx,1);
    else
        x0_pre = ad_est_acc(iSimStep-1,:)';
    end
    
    
    % b_ref generation
    if iSimStep == 1
        b_ref = zeros(N*nx,1);
    elseif iSimStep == 2
        b_ref = -M1*B*U_acc(iSimStep-1,:)';
    else
        b_ref = -M1*B*U_acc(iSimStep-1,:)' -M2*B*U_acc(iSimStep-2,:)';
    end
    
    
    % r, c generation
    r = 2*(B_conv)'*Q_tilda*(M1*x0 + M2*x0_pre + b_ref);
    c = (M1*x0 + M2*x0_pre + b_ref)'*Q_tilda*(M1*x0 + M2*x0_pre + b_ref); % constant
    
    %% MPC Simulate by CVX
    % Solving MPC via CVX
    tic
    cvx_precision best % cvx precision settings - low, medium, default, high, best options are exist.
    cvx_begin
    cvx_solver SDPT3 % sedumi, mosek
    variable U5(nu*N,1)
    subject to
    U5 <= U_max;
    U_min <= U5;
    E*U5 <= dU_max;
    dU_min <= E*U5;
    minimize( (U5)'*H*(U5) + r'*(U5) + c )
    cvx_end
    T_cvx(iSimStep,1) = toc;
    
    solver_objective(iSimStep,1) = cvx_optval;
    solver_feasibility{iSimStep,1} = cvx_status;

    %% Considering the input constraints
    U = U5;
    for i=1:size(U,1)
        if U(i,1) < 0
            U_v(i,1) = -( -coeff_b + sqrt((coeff_b)^2 - 4*coeff_a*U(i,1)*unit_change) )/(2*coeff_a);
        else
            U_v(i,1) = ( -coeff_b + sqrt((coeff_b)^2 + 4*coeff_a*U(i,1)*unit_change) )/(2*coeff_a);
        end
    end
    
    U_v_acc = [U_v_acc; U_v(1:nu)'];
    
    %% Phase correction by DM
    J = [J; U'*H*U + r'*U + c];
    u_prev = U(1:nu,1); % u_prev = U(nu*(N-1)+1:nu*N,1);
    ad_cor = B*u_prev; % ad_cor = A*x0 - A*ad_acc(iSimStep,:)' + B*u_prev;
    
    X_predicted = M1*x0 + M2*x0_pre + B_conv*U + b_ref;
    
    x_prev = X_predicted(1:nx,1);
    
    Zs_cor = zeros(nx,len,len);
    for j1=1:nx
        Zs_cor(j1,:,:) = ad_cor(j1,1).*squeeze(Zs(j1+1,:,:)); % piston element is removed
    end
    Zs_cor = squeeze(sum(Zs_cor,1));
    
    phase_cor(:,:,iSimStep) = Zs_cor;
    
    %%Iteration
    
    for i1=1:N
        X_err(iSimStep,i1) = norm((X_predicted((i1-1)*nx+1:i1*nx,1)),2);
    end
    
    X_err_low = [X_err_low; norm(x_prev(1:nx_est,1))];
    X_err_high = [X_err_high; norm(x_prev(nx_est+1:end,1))];
    
    
    if iSimStep == 1
        du_prev = u_prev - zeros(nu,1);
    else
        du_prev = u_prev - U_acc(iSimStep-1,:)';
    end
        
    ad_cor_acc = [ad_cor_acc; ad_cor'];
    U_acc = [U_acc; u_prev'];
    dU_acc = [dU_acc; du_prev'];
    
    X_acc = [X_acc; x_prev'];
    X_acc_err = [X_acc_err; norm(x_prev)];
    
    RMS_phase_res(iSimStep,1) = rms(rms(phase_res(:,:,iSimStep)));
    
    T_sim(iSimStep,1) = toc;
    disp(iSimStep)
    
end

%% System Matrix design (MPC_DesignMatrices_v4)
% Linear MPC의 state equation 수정
% Delta U constraint 생성에 필요한 matrx, E 추가
function [H, M1, M2, Q_tilda, R_tilda, B_conv, B_conv_pre1, B_conv_pre2, E] = MPC_DesignMatrices_v4(A1,A2,B,nx,nu,N,Q,P,R)

M1 = zeros(nx*N,nx); M2 = zeros(nx*N,nx);
for i=0:N-1
    if i==0
        M1(nx*i+1:nx*(i+1),:) = A1;
        M2(nx*i+1:nx*(i+1),:) = A2;
    elseif i==1
        M1(nx*i+1:nx*(i+1),:) = (A1)^2 + A2;
        M2(nx*i+1:nx*(i+1),:) = A1*A2;
    else
        M1(nx*i+1:nx*(i+1),:) = A1*(M1(nx*(i-1)+1:nx*i,:)) + A2*(M1(nx*(i-2)+1:nx*(i-1),:));
        M2(nx*i+1:nx*(i+1),:) = M1(nx*(i-1)+1:nx*i,:)*A2;
    end
end

B_cell = repmat({B}, 1, N);
B_conv = blkdiag(B_cell{:});

B_conv_pre1 = -M1*B;
B_conv_pre2 = -M2*B;


Q_tilda = zeros(nx*N,nx*N);
for i=0:N-2 
    Q_tilda((i*nx)+1:(i+1)*nx,(i*nx)+1:(i+1)*nx) = Q;
end
Q_tilda((N-1)*nx+1:end,(N-1)*nx+1:end) = P;


R_tilda = zeros(nu*N,nu*N);
for i=0:N-1
    R_tilda((i*nu)+1:(i+1)*nu,(i*nu)+1:(i+1)*nu) = R;
end


if N==1
    E = eye(nu,nu);
else
    E = zeros(nu*N,nu*N);
    for i=0:N-1
        if i==0
            E(i*nu+1:(i+1)*nu,i*nu+1:(i+1)*nu) = eye(nu,nu);
        else
            E(i*nu+1:(i+1)*nu,(i-1)*nu+1:i*nu) = -eye(nu,nu);
            E(i*nu+1:(i+1)*nu,i*nu+1:(i+1)*nu) = eye(nu,nu);
        end
    end
end

% matrix for cost function
H = 0.5*((B_conv)'*Q_tilda*B_conv + ((B_conv)'*Q_tilda*B_conv)') + R_tilda ; %J = U'*H*U + r'*U + c

end
