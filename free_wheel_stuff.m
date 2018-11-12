%[0m[ INFO] [1542024719.388820084]: Data in format: t, w_ref, v_ref, L_vel, R_vel, w_curent, v_current, Pl, Pr, al, ar[0m
%[0m[ INFO] [1542024719.389076959]: 0.050000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000[0m

[~, ~, ~, t,wr, vr, LV, RV, w, v, PL, PR, aL, aR, ~] = textread(sprintf('free_wheel_rand.txt'),'%s %s %s %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f %s');


figure();
hold on;
plot(t,vr);
plot(t,v);
legend('vr','v');
hold off;

figure();
hold on;
plot(t,LV);
plot(t,v);
legend('LV','v');
hold off;

