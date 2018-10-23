%% GEt data
[t5,LP5, RP5, LV5, RV5] = textread('50.txt','%f,%f,%f,%f,%f');
[t10,LP10, RP10, LV10, RV10] = textread('100.txt','%f,%f,%f,%f,%f');
[t2,LP2, RP2, LV2, RV2] = textread('20.txt','%f,%f,%f,%f,%f');

[t100,LP100, RP100, LV100, RV100] = textread('100olov.txt','%f,%f,%f,%f,%f');

%%
figure(50);
hold on;
plot(t5,LV5);
plot(t5,RV5);
plot(t5,LP5/100);
plot(t5,RP5/100);
tao50 = 0.5;
gain = 0.95;
[y,t] = step(tf([gain*tao50],[1,tao50]),20);
plot(t+3,0.5*y);
legend('LV','RV','LP','RP');
hold off;

%% 100 percent but too short
figure(10);
hold on;
plot(t10,LV10);
plot(t10,RV10);
plot(t10,LP10/100);
plot(t10,RP10/100);
legend('LV','RV','LP','RP');
hold off;


%% yo2
figure(20);
hold on;
plot(t2,LV2);
plot(t2,RV2);
plot(t2,LP2/100);
plot(t2,RP2/100);
tao20 = 0.5;
gain = 0.5
[y,t] = step(tf([gain*tao20],[1,tao20]),20);
plot(t+3,0.2*y);
legend('LV','RV','LP','RP');
hold off;

%% yo1
figure(100);
hold on;
plot(t100,LV100);
plot(t100,RV100);
plot(t100,LP100/100);
plot(t100,RP100/100);
tao100 = 0.5;
gain100 = 1.2;
[y100,to100] = step(tf([gain100*tao100],[1,tao100]),40);
plot(to100+3,y100);
legend('LV','RV','LP','RP');
hold off;

figure(101);
hold on;
plot(0.5.*(t100(1:end-1)+t100(2:end)),diff(LV100)/diff(t100));
plot(0.5.*(t100(1:end-1)+t100(2:end)),diff(RV100)/diff(t100));
plot(0.5.*(to100(1:end-1)+to100(2:end))+3,diff(y100)/diff(to100))
%plot(t5,diff(RV5));
%plot(t5,LP5/100);
%plot(t5,RP5/100);
legend('LV','RV','LP','RP');
hold off;

%%

figure(69);
hold on;
subplot(2,1,1);
semilogx(20*log10(abs(fft(LP100/100).\fft(LV100))));
subplot(2,1,2);
semilogx(phase(fft(LP100/100).\fft(LV100)));



hold off;