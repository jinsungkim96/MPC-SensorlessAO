function [P,h] = fast_mpc_ineq_const(obj)

%% Initial condition check
if (size(obj.x_min,1)~=size(obj.Q,1)) || (size(obj.x_max,1)~=size(obj.Q,1))
    error('Check the state inequality constraints dimensions');
end
if (size(obj.u_min,1)~=size(obj.R,1)) || (size(obj.u_max,1)~=size(obj.R,1))
    error('Check cotrol iequality constraint dimension');
end

%% Parameters
x_min = obj.x_min;
x_max = obj.x_max;
u_min = obj.u_min;
u_max = obj.u_max;
du_min = obj.du_min;
du_max = obj.du_max;
u_prev = obj.u_prev;
T = obj.T;
n = size(x_min,1);
m = size(u_min,1);
x = obj.x0;

%% Inequality constraint construction
% P = zeros(2*T*(n+m),T*(n+m));
% h = zeros(2*T*(n+m),1);

% for i=1:2*(m+n):size(P,1)-2*(m+n)+1
%     if i==1
%         P(i:i+2*(m+n)-1,i:i+(m+n)-1) = [eye(m) zeros(m,n);(-eye(m)) zeros(m,n);...
%             zeros(n,m) eye(n); zeros(n,m) -(eye(n))];
%     else
%         P(i:i+2*(m+n)-1,(i+1)/2:(i+1)/2+(m+n)-1) = [eye(m) zeros(m,n);(-eye(m)) zeros(m,n);...
%             zeros(n,m) eye(n); zeros(n,m) -(eye(n))];
%     end
% end

% for i=1:2*(m+n):2*T*(m+n)-2*(m+n)+1
%    h(i:2*(m+n)+i-1) = [u_max;-u_min;x_max;-x_min];
% end

P_box = zeros(2*T*(m),T*(n+m));
h_box = zeros(2*T*(m),1);

%% box constraint
for i=1:2*(m):2*T*(m)-2*(m)+1
    if i==1
        P_box(i:i+2*(m)-1,i:i+(m)-1) = [eye(m); -eye(m)];
    else
        P_box(i:i+2*(m)-1,(i+1)/2+n*((i-1)/(2*m)):(i+1)/2+(m+n*((i-1)/(2*m)))-1) = [eye(m); -eye(m)];
    end
end

for i=1:2*(m):2*T*(m)-2*(m)+1
   h_box(i:2*(m)+i-1) = [u_max;-u_min];
end

P_ramp = [];
h_ramp = [];

%% ramp constraint
% P_ramp = zeros(2*T*(m),T*(n+m));
% h_ramp = zeros(2*T*(m),1);
% 
% for i=1:2*(m):2*T*(m)-2*(m)+1
%     if i==1
%         P_ramp(i:i+2*(m)-1,i:i+(m)-1) = [eye(m); -eye(m)];
%     else
%         P_ramp(i:i+2*(m)-1,(m+n)*(((i-1)/(2*m))-1)+1:(m+n)*((i-1)/(2*m))+m) = [-eye(m) zeros(m,n) eye(m); eye(m) zeros(m,n) -eye(m)];
%     end
% end
% 
% for i=1:2*(m):2*T*(m)-2*(m)+1
%     if i==1
%         h_ramp(i:2*(m)+i-1) = [u_prev+du_max; -u_prev-du_min];
%     else
%         h_ramp(i:2*(m)+i-1) = [du_max; -du_min];
%     end
% end

P = [P_box; P_ramp];
h = [h_box; h_ramp];

end