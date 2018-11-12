%[0m[ INFO] [1542024719.388820084]: Data in format: t, w_ref, v_ref, L_vel, R_vel, w_curent, v_current, Pl, Pr, al, ar[0m
%[0m[ INFO] [1542024719.389076959]: 0.050000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
% 0.000000, 0.000000, 0.000000, 0.000000[0m

[~, ~, ~, t,wr, vr, LV, RV, w, v, PL, PR, aL, aR, ~] = textread(sprintf('free_wheel_rand.txt'), ...
    '%s %s %s %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f %s');

testInt = 1:2000;
valInt = 2001:3000;
LVtest = LV(testInt);

%% try arx

%T = T- V- P-

% guess params
%xhat = [0.8, -0.5, 1];
xhat = [0.020687152912531  -0.100178819223951   0.108139219408348]; %from many runs
y = normrnd(zeros(numel(t(testInt)),1),1); %zeros(numel(t),1);

%y = [0;diff(LV(testInt))./diff(t(testInt))]; % guess that tao equal to accelleration.
phi = PL(testInt)./100; % Normalized

%u = [y, LV(1:end-1), PL(1:end-1)]; % y doesn't need offset here since all zeros anyways.

%y = xhat*u;
% estimate params 
%xhat = y/u;
figure();
hold on; plot([0;diff(LVtest(testInt))./diff(t(testInt))]); hold off;

%% loop!
m = 200;
options = optimoptions('fmincon','Display','off');
for i=1:m

    fx = @(x)taomin(y,x,LVtest,phi); %function handle for funciton to minimize
    [xhat,~] = fmincon(fx,xhat,[],[],[],[],[0,-inf,0],[1,0,inf],[],options); % Minimize!
    
    fy = @(y)taomin(y,xhat,LVtest,phi); %function handle for funciton to minimize
    [y,~] = fmincon(fy,y,[],[],[],[],[],[],[],options); % Minimize!
    %hold on;
    %plot(y);
    %hold off;
    %drawnow;
    i
    xhat

end

%% Validate

figure();
hold on;
plot([0;diff(LV(testInt))./diff(t(testInt))],'k');
plot(xhat*[y,LV(testInt), PL(testInt)]','r');
plot(2*pi*LV(testInt),'c');
plot(PL(testInt)/5,'m');

legend('acc','estimated tao','LV','PL');

hold off;


