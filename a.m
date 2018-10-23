

%[0m[ INFO] [1539769436.604612613]: 0.010000,0.000000,0.000000,0.000000,0.000000[0m
%%
pwrs = [15:5:100];
t_s = 200;
t_e = 200+400;
t_end = (2+8)*100;
N = numel(pwrs);
gains = zeros(N,2); % R first then L
taos = zeros(N,2);
m = 103.5;
r = 0.165;

V_max = 13;
V_steps = V_max*(pwrs/100);

Ri = 1.3; % Internal R of engines, same for both.
Ki_R = Ri/0.4308;
Ki_L = Ri/0.4639;

Kws = zeros(N,2);
tunes = zeros(N,2);

%%

maxacc = zeros(N,2);

for i=1:N
    pwr = pwrs(i);
    V = V_steps(i);
    [~, ~, ~, t,LP, RP, LV, RV, ~] = textread(sprintf('%dpwr.txt',pwr),'%s %s %s %f,%f,%f,%f,%f %s');
    t = t(1:t_e-t_s+1);
    LP = LP(t_s:t_e);
    RP = RP(t_s:t_e);
    LV = LV(t_s:t_e);
    RV = RV(t_s:t_e);
    
    
    x0 = [1.1,2];             %Initial paramerer guesses.
    flhR = @(x) stepmin(x,RV,t,pwr);        %function handle for funciton to minimize
    flhL = @(x) stepmin(x,LV,t,pwr);        %function handle for funciton to minimize
    [xR,~] = fmincon(flhR,x0,[],[],[],[],[0,0],[inf,inf]); % Minimize!
    [xL,~] = fmincon(flhL,x0,[],[],[],[],[0,0],[inf,inf]); % Minimize!
    gainR = xR(1);
    taoR = xR(2);
    gainL = xL(1);
    taoL = xL(2);
    
    gains(i,:) = [gainR,gainL];
    taos(i,:) = [taoR,taoL];
    
    sysR = tf(gainR/taoR,[1,1/taoR]);
    sysL = tf(gainL/taoL,[1,1/taoL]);
    opt = stepDataOptions('StepAmplitude',pwr/100);
    [yR,~] = step( sysR ,t,opt);
    [yL,~] = step( sysL ,t,opt);
    
    f = figure('rend','painters','pos',[910 610 900 600]);
    hold on;
    title(sprintf('%d pwr',pwr));
    plot(t,LV);
    plot(t,RV);
    plot(t,yR);
    plot(t,yL);
    hold off;
    close(f);
    
    accR = diff(yR)./diff(t);
    accL = diff(yL)./diff(t);
    acc_t = (t(2:end)+t(1:end-1))/2;
    
    f = figure('rend','painters','pos',[10 10 900 600]);
    hold on;
    plot( acc_t , accR); % R
    plot( acc_t , accL); % L
    maxacc(i,:) = [max(accR), max(accL)];
    legend('R','L');
    title('Accelleration curve');
    hold off;
    close(f);
    
    
    %tuneL = 2.2;
    %tuneR = 1.8;
    %Kw_L = (r*V - tuneL*r^2*0.5*accL.*Ri*m/Ki_L)./( 0.5*(yL(1:end-1)+yL(2:end)) ); % Half of mass due to two motors helping accelleration
    %Kw_R = (r*V - tuneR*r^2*0.5*accR.*Ri*m/Ki_R)./( 0.5*(yR(1:end-1)+yR(2:end)) );
    
    Kw_L_data = @(tune) (r*V - tune*r^2*accL.*Ri*m/Ki_L)./( 0.5*(yL(1:end-1)+yL(2:end)));
    Kw_R_data = @(tune) (r*V - tune*r^2*accR.*Ri*m/Ki_R)./( 0.5*(yR(1:end-1)+yR(2:end)));
    
    fKwL = @(tune) mean((Kw_L_data(tune) - mean(Kw_L_data(tune))).^2);
    fKwR = @(tune) mean((Kw_R_data(tune) - mean(Kw_R_data(tune))).^2);
    
    [tune_R, fvalR] = fmincon(fKwR, 2);
    [tune_L, fvalL] = fmincon(fKwL, 2);
    
    Kw_R = mean(Kw_R_data(tune_R));
    Kw_L = mean(Kw_L_data(tune_L));
    Kws(i,:) = [Kw_R, Kw_L];
    tunes(i,:) = [tune_R, tune_L];
    
    f = figure('rend','painters','pos',[200 200 900 600]);
    hold on;
    offset = 50;
    KWRD = Kw_R_data(tune_R);
    KWLD = Kw_L_data(tune_L);
    
    plot( acc_t(offset:end) , KWRD(offset:end)); % R
    plot( acc_t(offset:end) , KWLD(offset:end)); % L
    legend('R','L');
    title(sprintf('Kw curve for step power: %d',pwr));
    hold off;
    %close(f);
    
end
%close all;
figure();
plot(Kws);
title('Kws');
figure();
plot(tunes);
title('tunes');
%%

figure();
hold on;
scatter(pwrs,gains(:,2),300,'x'); % L
scatter(pwrs,gains(:,1),300,'.'); % R
legend('L gains','R gains');
hold off;
figure();
hold on;
scatter(pwrs(1:end),taos(1:end,2),300,'x') % L
scatter(pwrs(1:end),taos(1:end,1),300,'.') % R
title('taos');
legend('taos L','taos R');
hold off;


figure();
hold on;
scatter(pwrs,maxacc(:,2),300,'x') % L
scatter(pwrs,maxacc(:,1),300,'.') % R
title('maxacc');
legend('L','R');
hold off;

figure();
hold on;
scatter(pwrs,maxacc(:,2).*103.5*0.165,300,'x') % L
scatter(pwrs,maxacc(:,1).*103.5*0.165,300,'.') % R
title('max torque');
legend('L','R');
hold off;

%% Plot U vs max torque 

UvsTorqueL = V_steps'./(maxacc(:,2).*m*r);
UvsTorqueR = V_steps'./(maxacc(:,1).*m*r);

figure();
hold on;
scatter(pwrs,UvsTorqueL,300,'x') % L
scatter(pwrs,UvsTorqueR,300,'.') % R
title('U/max torque = Ri/Ki');
legend('L','R');
hold off;



%%
pwr = 75;
t_s = 200;
t_e = 200+800;
[~, ~, ~, t,LP, RP, LV, RV, ~] = textread(sprintf('%dpwr.txt',pwr),'%s %s %s %f,%f,%f,%f,%f %s');

t = t(1:t_e-t_s+1);
LP = LP(t_s:t_e);
RP = RP(t_s:t_e);
LV = LV(t_s:t_e);
RV = RV(t_s:t_e);


x0 = [1.1,2];             %Initial paramerer guesses.
flh = @(x) stepmin(x,RV,t,pwr);        %function handle for funciton to minimize
[x,fval] = fmincon(flh,x0,[],[],[],[],[0,0],[inf,inf]); % Minimize!
fval    %Resulting function value
x       %Parameters

%%
gain = x(1); %1.1;
tao = x(2); %2;
sys = tf(gain*tao,[1,tao]);
opt = stepDataOptions('StepAmplitude',pwr/100);
[y,~] = step( sys ,t,opt);

figure('rend','painters','pos',[10 10 900 600])
hold on;
title(sprintf('%d pwr',pwr));
plot(t,LV)
plot(t,RV);
plot(t,y);

hold off;