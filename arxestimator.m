% Assumes youve run a.m and then efteanal.m

u;
bsize = 2^6;
LVs_smooth = conv(LVs,bartlett(bsize)./sum(bartlett(bsize)),'same');


acc_LVss = diff(LVs_smooth)./diff(ts);
t_acc = 0.5*(ts(1:end-1) + ts(2:end));
ws = 0.5*(LVs(1:end-1) + LVs(2:end));
us = V_max*0.5*(u(1:end-1)+u(2:end));

%plot(0.5*(ts(1:end-1) + ts(2:end)),diff(yf)./diff(ts));

%% ARX 1,1,1
B1 = acc_LVss(2:end)'/[us(1:end-1)'; ws(1:end-1)'];

figure();
hold on;

plot(t_acc,acc_LVss);
plot(t_acc(2:end), B1*[us(1:end-1)'; ws(1:end-1)'])
legend('actual data','estimated acc');

hold off;

%% ARX 2,1,2

Input = [acc_LVss(2:end-1)'; us(2:end-1)'; ws(2:end-1)'; ws(1:end-2)'];
B2 = acc_LVss(3:end)'/Input;

figure();
hold on;

plot(t_acc,acc_LVss);
plot(t_acc(3:end), B2*Input)
legend('actual data','estimated acc');

hold off;

%%


Input = [acc_LVss(1:end-1)'; us(1:end-1)'; ws(1:end-1)'];
B3 = acc_LVss(2:end)'/Input

figure();
hold on;

plot(t_acc,acc_LVss);
plot(t_acc(2:end), B3*Input)
legend('actual data','estimated acc');

hold off;
%%


Input = [acc_LVss(1:end-2)'; acc_LVss(2:end-1)'; us(2:end-1)'; ws(2:end-1)'];
B4 = acc_LVss(3:end)'/Input

figure();
hold on;

plot(t_acc,acc_LVss);
plot(t_acc(3:end), B4*Input)
legend('actual data','estimated acc');

hold off;

%%
figure();
hold on;

plot(t_acc,ws);
plot(t_acc(3:end), 0.01*cumsum(B4*Input))
legend('actual data','estimated acc');

hold off;

%% Try slow sampling
nskip = 5;
Input = [acc_LVss(1:nskip:end-2)'; acc_LVss(2:nskip:end-1)'; us(2:nskip:end-1)'; ws(2:nskip:end-1)'];
B4 = acc_LVss(3:nskip:end)'/Input

figure();
hold on;

plot(t_acc,acc_LVss);
plot(t_acc(3:nskip:end), B4*Input)
legend('actual data','estimated acc');

hold off;

%%
figure();
hold on;

plot(t_acc,ws);
plot(t_acc(3:nskip:end), 0.05*cumsum(B4*Input))
legend('actual data','estimated acc');

hold off;
