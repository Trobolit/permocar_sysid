%%
[~, ~, ~, t,LP, RP, LV, RV, ~] = textread(sprintf('megadata.txt'),'%s %s %s %f,%f,%f,%f,%f %s');

V_max = 12.8;

LV(270:340) = 0;
RV(270:340) = 0;

LV(535:620) = 0;
RV(535:620) = 0;

LV(817:924) = 0;
RV(817:924) = 0;

LV(1231:1394) = 0;
RV(1231:1394) = 0;

LV(1609:1761) = 0;
RV(1609:1761) = 0;

LV(2098:2571) = 0;
RV(2098:2571) = 0;
LP(2098:2571) = 0;
RP(2098:2571) = 0;

LV(2713:2817) = 0;
RV(2713:2817) = 0;

LV(3367:3571) = 0;
RV(3367:3571) = 0;

LV(4101:4214) = 0;
RV(4101:4214) = 0;

LV(4433:4509) = 0;
RV(4433:4509) = 0;

LV(4674:4803) = 0;
RV(4674:4803) = 0;

LV(4962:end) = 0;
RV(4962:end) = 0;