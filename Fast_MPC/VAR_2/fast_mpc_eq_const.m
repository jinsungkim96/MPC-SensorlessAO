function [C,b] = fast_mpc_eq_const(obj)

%% Parameters
A1 = obj.A1;
A2 = obj.A2;
B = obj.B;
w = obj.w;
n = size(A1,2);
m = size(B,2);
x0 = obj.x0;
x0_pre = obj.x0_pre;
u = obj.R;
T = obj.T;
C = zeros(T*n,T*(n+m));
b = zeros(T*n,1);
xf = obj.x_final;

%% Initial condition check
if isempty(A1)
    error('Define the state dynamics/equality constrained matrix');
elseif isempty(A2)
    error('Define the state dynamics/equality constrained matrix');
elseif isempty(B)
    error('Define the control dynamics/equality constrained matrix');
end

if size(A1,2) ~= size(x0,1)
    error('The equality state dynamics matrix size does not match');
elseif size(A2,2) ~= size(x0_pre,1)
    error('The equality state dynamics matrix size does not match');
elseif size(B,2) ~= size(u,2)
    error('The equality control dynamics matrix size does not match');
elseif isempty(w)
    w = zeros(n,1);
end

%% Equality constraint construction
C(1:n,1:m+n) = [-B eye(n)];
b(1:n) = A1*x0 + A2*x0_pre + w(1:n);

for i=1:T-1
    if i==1
        C(n*i+1:n*(i+1),m+1:(n+m)*(i+1)) = [-A1 -B eye(n)];
        b(n*i+1:n*(i+1)) = A2*x0 + w(n*i+1:n*(i+1));
    else
        C(n*i+1:n*(i+1),m+1+(n+m)*(i-2):(n+m)*(i+1)) = [-A2 zeros(n,m) -A1 -B eye(n)];
        b(n*i+1:n*(i+1)) = w(n*i+1:n*(i+1));
    end
end


% for i=n:n:T*n-n+1
%     if i==n
%         C(i+1:i+n,i+((i/n)-1)*(n+m+n):i+n+((i/n)-1)*(n+m+n)+m+n-1) = [-A1 -B eye(n)];
%         b(i+1:i+n) = A2*x0 + w(i+1:i+n);
%         % b(i+1:i+n) = w;
%     else
%         C(i+1:i+n,((i/n)-1)*(n+m)+m+1:((i/n)-1)*(n+m)+m+m+n+n) = [-A2 zeros(m) -A1 -B eye(n)];
%         b(i+1:i+n) = w(i+1:i+n);
%         % b(i+1:i+n) = w;
%     end
% end
% C(1:n,1:m+n) = [-B eye(n)];
% b(1:n) = A1*x0 + A2*x0_pre + w(1:n);

% b(end-n+1:end) = xf;
if isempty(xf)~=1
   b = [b;xf];
   C = [C;zeros(n,size(C,2))];
end
C(end-n+1:end,end-n+1:end) = eye(n);
end