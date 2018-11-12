%[0m[ INFO] [1542024719.388820084]: Data in format: t, w_ref, v_ref, L_vel, R_vel, w_curent, v_current, Pl, Pr, al, ar[0m
%[0m[ INFO] [1542024719.389076959]: 0.050000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
% 0.000000, 0.000000, 0.000000, 0.000000[0m

[~, ~, ~, t,wr, vr, LV, RV, w, v, PL, PR, aL, aR, ~] = textread(sprintf('free_wheel_rand.txt'), ...
    '%s %s %s %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f %s');


%% try arx

%T = T- V- P-

% guess params
xhat = [0.8, -0.5, 1];
y = zeros(numel(t),1);

phi = PL./100; % Normalized

%u = [y, LV(1:end-1), PL(1:end-1)]; % y doesn't need offset here since all zeros anyways.

%y = xhat*u;
% estimate params 
%xhat = y/u;

m = 100;
for i=1:m

    fx = @(x)taomin(y,x,LV,phi); %function handle for funciton to minimize
    [xhat,~] = fmincon(fx,xhat,[],[],[],[],[0,0,0],[1,inf,inf]); % Minimize!
    
    fy = @(y)taomin(y,xhat,LV,phi); %function handle for funciton to minimize
    [y,~] = fmincon(fy,y,[],[],[],[],[],[]); % Minimize!
    
    i
    xhat

end
