%%
[~, ~, ~, t,LP, RP, LV, RV, ~] = textread(sprintf('megadata.txt'),'%s %s %s %f,%f,%f,%f,%f %s');

V_max = 12.8;

LV(270:340) = 0;
RV(270:340) = 0;

LV(535:620) = 0;
RV(535:620) = 0;

LV(817:924) = 0;
RV(817:924) = 0;

LV(1231:1394) = 0;
RV(1231:1394) = 0;

LV(1609:1761) = 0;
RV(1609:1761) = 0;

LV(2098:2571) = 0;
RV(2098:2571) = 0;
LP(2098:2571) = 0;
RP(2098:2571) = 0;

LV(2713:2817) = 0;
RV(2713:2817) = 0;

LV(3367:3571) = 0;
RV(3367:3571) = 0;

LV(4101:4214) = 0;
RV(4101:4214) = 0;

LV(4433:4509) = 0;
RV(4433:4509) = 0;

LV(4674:4803) = 0;
RV(4674:4803) = 0;

LV(4962:end) = 0;
RV(4962:end) = 0;

figure(1);
hold on;
plot(LP/100)
plot(LV);
plot(RV);
legend('LP','LV','RV');
hold off;


%%
u = LP;

bsize = 30; %40;
LV_smooth = conv(LV,bartlett(bsize)./sum(bartlett(bsize)),'same');
%LV_smooth = LV;
acc_LVs = diff(LV_smooth)./diff(t(1:end));
acc_LV = diff(LV)./diff(t(1:end));
t_acc = 0.5*(t(1:end-1) + t(2:end));
wL_acc = 0.5*(LV(1:end-1) + LV(2:end));
uL_acc = V_max*0.5*(u(1:end-1)+u(2:end));

Input = [acc_LVs(1:end-2)'; acc_LVs(2:end-1)'; uL_acc(2:end-1)'; wL_acc(2:end-1)'];
B = acc_LVs(3:end)'/Input;

figure(2);
hold on;


plot(t_acc,acc_LV);
plot(t_acc(3:end), B*Input)
legend('actual data','estimated acc');

hold off;


figure(3);
hold on;

plot(t_acc,wL_acc);
plot(t_acc(3:end), 0.05*cumsum(B*Input))
legend('actual data','estimated LV');

hold off;


%% ETFE
y = LV;
Ghathat = fft(u/100).\fft(y); % Straight forward.
figure();
hold on;
subplot(2,1,1);
semilogx(20*log10(abs(Ghathat))); %Custom bode plot
subplot(2,1,2);
semilogx(phase(Ghathat));
hold off;

%%

bsize = 2^6; % Assuming size of 2^n is the most efficient for underlying algorithms.
%Gsmooth = conv(bartlett(bsize),Ghathat) / sum(bartlett(bsize)); % has bad
%effects at endpoints of data due to zeros being averiaged together with
%small amounts of data. Below solves that.
Gsmooth = conv(bartlett(bsize),Ghathat.*abs(fft(u/100)).^2) ./ conv(bartlett(bsize),abs(fft(u/100)).^2);
figure();
hold on;
subplot(2,1,1);
semilogx(20*log10(abs(Gsmooth)));
subplot(2,1,2);
semilogx(phase(Gsmooth));
hold off;

%%

%[acc_LVs(1:end-2)'; acc_LVs(2:end-1)'; uL_acc(2:end-1)'; wL_acc(2:end-1)'];
Phi = (1/B(3)) * (acc_LVs(3:end) - B(1)*acc_LVs(1:end-2) - B(2)*acc_LVs(2:end-1) -B(4)*wL_acc(2:end-1));

figure();
hold on;

plot(t_acc(3:end),Phi/20);
plot(t,LP);

hold off;

%% estimate phi instead

LPa = 0.5*(LP(1:end-1) + LP(2:end));

Input1 = [acc_LV(1:end-2)'; acc_LV(2:end-1)'; acc_LV(3:end)'; wL_acc(3:end)'];
B1 = LPa(3:end)'/Input1;
Input2 = [acc_LV(1:end-1)'; acc_LV(2:end)'; wL_acc(2:end)'];
B2 = LPa(2:end)'/Input2;
Input3 = [acc_LV(1:end-1)'; acc_LV(2:end)'; wL_acc(2:end)'; LPa(1:end-1)'];
B3 = LPa(2:end)'/Input3;
%B3 = [1.15, 5.19, 5.08, 0.92]; % OVERRIDE (These are numerically prop what will be implemented, they work.)
Input4 = [acc_LV(2:end-1)'; acc_LV(3:end)'; wL_acc(3:end)'; LPa(2:end-1)'; LPa(1:end-2)'];
B4 = LPa(3:end)'/Input4;


figure();
hold on;


plot(t_acc,LPa);
plot(t_acc(3:end)', B1*Input1);
plot(t_acc(2:end)', B2*Input2);
plot(t_acc(2:end)', B3*Input3);
plot(t_acc(3:end)', B4*Input4);
legend('actual data','estimated acc extra','estimated simple','with LP mem', 'with phi mem');

hold off;


%% Estimate phi using data in a way we will have it

%{
acc_LV = diff(LV)./diff(t(1:end));
t_acc = 0.5*(t(1:end-1) + t(2:end));
wL_acc = 0.5*(LV(1:end-1) + LV(2:end));
uL_acc = V_max*0.5*(u(1:end-1)+u(2:end));
LPa = 0.5*(LP(1:end-1) + LP(2:end));
%}

% Input = [a_wanted, a_last=(v_now - v_last)/2, V_now=(2*v+a_wanted*dt)/2, Phi_last]
% Input = [a_wanted, v_now, v_last, Phi_last]
acc_LV = diff([0;LV])./diff( [t(1)-(t(2)-t(1));t] );
acc_LV = acc_LV(2:end); % remove a0 since it does not exist in reality.

% a_wanted, v, v-, LP-
Input5 = [acc_LV'; LV(1:end-1)'; [0;LV(1:end-2)]'; [0;LP(1:end-2)]' ];
B5 = LP(1:end-1)'/Input5;
% a_wanted, v, v-, LP-, a_wanted-
Input6 = [acc_LV'; LV(1:end-1)'; [0;LV(1:end-2)]'; [0;LP(1:end-2)]'; [0;acc_LV(1:end-1)]' ];
B6 = LP(1:end-1)'/Input6;

% a_wanted v, v-, P-
Input7 = [[acc_LV(1:end)]'; [LV(1:end-1)]'; [0;LV(1:end-2)]'; [0;LP(1:end-2)]' ];
B7 = LP(2:end)'/Input7;

% a_wanted v, P-
Input8 = [acc_LV'; LV(1:end-1)'; [0;LP(1:end-2)]' ];
B8 = LP(1:end-1)'/Input8;



figure();
hold on;

plot(t,LP);
plot(t(1:end-1)', B5*Input5);
%plot(t(1:end-1)', B6*Input6);
plot(t(1:end-1)', B7*Input7);
%plot(t_acc(2:end)', B3*Input3);
plot(t(1:end-1)', B8*Input8);
legend('actual data','estimated LP','in 7','simple AF');

hold off;

%%
figure();
hold on;

plot(t,LP);
plot(t(1:end-1)', B8*Input8);
y = nan(numel(t)-1,1);
for k=1:numel(t)-1
    if k>1
        y(k) = getPhi(acc_LV(k), LV(k), LP(k-1));
    else
        y(k) = getPhi(acc_LV(k), LV(k), 0 );
    end
end
%plot(t(1:end-1)', getPhi(acc_LV, LV(1:end-1), [0;LP(1:end-2)]));
plot(t(1:end-1)', y);
legend('actual data','simple AF','for loop');

hold off;
