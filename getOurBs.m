filterData;

acc_LV = diff([0;LV])./diff( [t(1)-(t(2)-t(1));t] );
acc_LV = acc_LV(2:end); % remove a0 since it does not exist in reality.
acc_RV = diff([0;RV])./diff( [t(1)-(t(2)-t(1));t] );
acc_RV = acc_RV(2:end); % remove a0 since it does not exist in reality.


%Left
% a_wanted v, P-
InputL = [acc_LV'; LV(1:end-1)'; [0;LP(1:end-2)]' ];
BL = LP(1:end-1)'/InputL

%Left
% a_wanted v, P-
InputR = [acc_RV'; RV(1:end-1)'; [0;RP(1:end-2)]' ];
BR = RP(1:end-1)'/InputR


figure();
hold on;

plot(t,LP);
plot(t(1:end-1)', BL*InputL);
plot(t(1:end-1)', BR*InputR);
legend('actual data','estimated LP','estimated RP');

hold off;
