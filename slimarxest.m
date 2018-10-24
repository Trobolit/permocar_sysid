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
