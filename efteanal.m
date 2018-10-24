% assumes a.m has been run.

% Stitch together all stepinfo data

n = 41000;

LPs = nan(n,1);
RPs = nan(n,1);
LVs = nan(n,1);
RVs = nan(n,1);
ts = nan(n,1);
m = 1;

for i=1:N
    pwr = pwrs(i);
    V = V_steps(i);
    [~, ~, ~, t,LP, RP, LV, RV, ~] = textread(sprintf('%dpwr.txt',pwr),'%s %s %s %f,%f,%f,%f,%f %s');
    %fprintf('sizes :: LP:%d, RP %d, LV %d, RV %d, t %d\n',size(LP),size(RP),size(LV),size(RV),size(t));
    p = numel(t);
    LPs(m:m+p-1,1) = LP;
    RPs(m:m+p-1,1) = RP;
    LVs(m:m+p-1,1) = LV;
    RVs(m:m+p-1,1) = RV;
    if m==1
        ts(m:m+p-1) = t;
    else
        ts(m:m+p-1) = t+ts(m-2);
    end
    m = m +p;
    
end

%%
figure();
hold on;

plot(ts,LPs/100);
plot(ts,LVs);

hold off;

%% ETFE
y = LVs;
u = LPs/100;
Ghathat = fft(u).\fft(y); % Straight forward.
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
Gsmooth = conv(bartlett(bsize),Ghathat.*abs(fft(u)).^2) ./ conv(bartlett(bsize),abs(fft(u)).^2);
figure();
hold on;
subplot(2,1,1);
semilogx(20*log10(abs(Gsmooth)));
subplot(2,1,2);
semilogx(phase(Gsmooth));
hold off;

