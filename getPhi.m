function [Phi] = getPhi(a_wanted, v_now, Phi_last)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

B = [9.611797484276181  11.449634200880601   0.846745613948956];

%a = tao/(m*r)


%size(a_wanted)
%size(v_now)
%size(Phi_last)
% a_wanted v, P-
%Phi = B*[a_wanted'; v_now'; Phi_last' ];

%{
if abs(Phi_last) > 100
    Phi_last = 100*sign(Phi_last);
end
%}

Phi = 9.61*a_wanted' + 11.45*v_now + 0.85*Phi_last;

if abs(Phi) > 100
    Phi = 100*sign(Phi);
end

end

