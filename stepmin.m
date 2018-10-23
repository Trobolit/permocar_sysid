function mymse = stepmin(x,y,t,pwr)
%STEPMIN Summary of this function goes here
%   Detailed explanation goes here
sys = tf(x(1)*1/x(2),[1,1/x(2)]);
opt = stepDataOptions('StepAmplitude',pwr/100);
[y_step,~] = step( sys ,t,opt);

mymse = mean((y-y_step).^2);

end

