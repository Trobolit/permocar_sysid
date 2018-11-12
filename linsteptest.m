%NFO] [1541420022.935877023]: Data in format: t, w_ref, v_ref, L_vel, R_vel
[~, ~, ~, t,wr, vr, LV, RV, ~] = textread(sprintf('linsteptest1.txt'),'%s %s %s %f,%f,%f,%f,%f %s');

%NFO] [1541424576.931526010]: Data in format: t, w_ref, v_ref, L_vel, R_vel, w_curent, v_current, Pl, Pr, al, ar
[~, ~, ~, t2,wr2, vr2, LV2, RV2, w2, v2, PL2, PR2, aL, aR, ~] = textread(sprintf('linsteptest2.txt'),'%s %s %s %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f %s');

%%

sys = tf([1.3],[1, 1.3]);
opt = stepDataOptions('StepAmplitude',0.82);
[y, T] = step(sys,6, opt);

figure();
hold on;
plot(T+1.4, y);
plot(t,-vr/2);
plot(t,LV);
plot(t,RV);
legend('sim','vr','LV','RV');
hold off;

%%

sys = tf([1.5],[1, 1.5]);
opt = stepDataOptions('StepAmplitude',0.82);
[y, T] = step(sys,6, opt);

figure();
hold on;
plot(T+1.5, y);
plot(t,-vr/2);
plot(t,LV);
plot(t,RV);
legend('sim','vr','LV','RV');
hold off;

%% 2
sys = tf([1.5],[1, 1.5]);
opt = stepDataOptions('StepAmplitude',0.82);
[y, T] = step(sys,6, opt);

figure();
hold on;
plot(T+4.9, y);
plot(t2,vr2);
plot(t2,v2);
plot(t2, aL);
plot(t2, PL2/100);
legend('sim','vr','v', 'aL_wanted', 'PL');
hold off;