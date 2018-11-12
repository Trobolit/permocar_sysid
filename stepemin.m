function [e] = stepemin(x, t, v)
%STEPEMIN Summary of this function goes here
%   Detailed explanation goes here

t = t(16:200);
v = v(16:200);
%x

ka = x(1) +x(3)*1i;
kb = x(2) -x(3)*1i;
kc = x(4);
% (s+p1)(s+p2)(s+p3) = (s+p1)(s^2 +(p2+p3)s + p2*p3) =
% s^3 +(p2+p3)s^2 + p2*p3*s + p1*s^2 +p1*(p2+p3)s + p1*p2*p3
% = s^3 + (p1+p2+p3)s^2 + p2*p34*s +p1*p2*s +p1*p3*s + p1*p2*p3
sys3 = tf([ka*kb*kc],[1, ka+kb+kc, ka*kb+ka*kc+kb*kc, ka*kb*kc]);
opt = stepDataOptions('StepAmplitude',x(5));
[y3, ~] = step(sys3,t-t(1), opt);

%plot(T3+3.35, y3); %3.5x
%plot(t,vr);

e = abs(sum((y3-v).^2));
%e
end

