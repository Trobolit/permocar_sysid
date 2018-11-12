
%NFO] [1541426926.251812080]: Data in format: t, w_ref, v_ref, L_vel, R_vel, w_curent, v_current, Pl, Pr, al, ar
[~, ~, ~, t,wr, vr, LV, RV, w, v, PL, PR, aL, aR, ~] = textread(sprintf('k1p5_v5.txt'),'%s %s %s %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f %s');

%%
x0 = [1.5, 1.5, 1.2, 1.5, 0.9];
fm = @(x)stepemin(x,t,v);        %function handle for funciton to minimize
[xhat,fval] = fmincon(fm,x0,[],[],[],[],[0,0,0,0,0],[inf,inf,inf,inf,1]); % Minimize!

%%

% (s+p1)(s+p2) = s^2 + p1+p2 + p1*p2
k1 = 1.5 + 1.1i; % 1.5
k2 = 1.5 - 1.1i;
sys = tf([k1*k2],[1, k1+k2, k1*k2]);



%ka = 1.5+1.2i;
%kb = 1.5-1.2i;
%kc = 1.5;
% (s+p1)(s+p2)(s+p3) = (s+p1)(s^2 +(p2+p3)s + p2*p3) =
% s^3 +(p2+p3)s^2 + p2*p3*s + p1*s^2 +p1*(p2+p3)s + p1*p2*p3
% = s^3 + (p1+p2+p3)s^2 + p2*p34*s +p1*p2*s +p1*p3*s + p1*p2*p3

ka = xhat(1) + xhat(3)*1i;
kb = xhat(2) - xhat(3)*1i;
bc = xhat(4);
sys3 = tf([ka*kb*kc],[1, ka+kb+kc, ka*kb+ka*kc+kb*kc, ka*kb*kc]);

opt = stepDataOptions('StepAmplitude',xhat(5));
[y3, T3] = step(sys,t-t(1), opt);

%%
figure();
hold on;
plot(T3+3.35, y3); %3.5x
plot(t,vr);
plot(t,v);
plot(t, aL);
plot(t, PL/100);
xlim([3,10]);
ylim([-0.1,1.2]);
legend('sim','vr','v', 'aL_wanted', 'PL');
hold off;