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
figure(456);
hold on;

plot(ts,LPs/100);
plot(ts,LVs);
bsize = 2^7;
LVs_filt = conv(LVs,bartlett(bsize)./sum(bartlett(bsize)),'same');
plot(ts,LVs_filt);

hold off;

%% ETFE

y = LVs;%_filt;
u = LPs/100;

Fs= 100; % insert here your frequency sampling in Hz
L=length(u); 
NFFT = 2^nextpow2(L);
f = Fs/2*linspace(0,1,NFFT/2+1);
%plot(f,2*abs(Y(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')

Ghathat = (fft(u,NFFT)/L).\(fft(y,NFFT)/L); % Straight forward.
figure();
hold on;
subplot(2,1,1);
semilogx(20*log10(abs(Ghathat(1:NFFT/2+1)))); %Custom bode plot
subplot(2,1,2);
semilogx(phase(Ghathat(1:NFFT/2+1)));
hold off;

%%

%Ghathat = Ghathat(1:NFFT/2+1);

bsize = 2^7; % Assuming size of 2^n is the most efficient for underlying algorithms.
%Gsmooth = conv(bartlett(bsize),Ghathat) / sum(bartlett(bsize)); % has bad
%effects at endpoints of data due to zeros being averiaged together with
%small amounts of data. Below solves that.
Gsmooth = conv(bartlett(bsize),Ghathat.*abs(fft(u,NFFT)/L).^2) ./ conv(bartlett(bsize),abs(fft(u,NFFT)/L).^2);
figure();
hold on;
subplot(2,1,1);
semilogx(20*log10(abs(Gsmooth(1:NFFT/2+1))));
subplot(2,1,2);
semilogx(phase(Gsmooth(1:NFFT/2+1)));
hold off;
