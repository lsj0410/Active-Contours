% Code only works for Nc = 4 (bc. Hadamard decoding process)
% Minor issue with findpeaks

clc; clear all;
Nc = 4;
z = zeros(Nc, 4, 128, 128);
c = 3*10^8;
f_c = 3*10^9; lambda = c/f_c;
B = 1.5*10^13; K = 3*10^15; f_s = 10^9;
T_c = B/K; T_0 = T_c*size(z, 1);
rmax = T_c*c/2/B*f_s; vmax = lambda/4/T_0;
dr = lambda/2; dt = dr*size(z, 2);
cnt = 1024;

% generating random input
% R = rand(3, 1)*rmax;
% v = (rand(3, 1) - 0.5)*vmax/Nc;
% a = (rand(3, 1) - 0.5)*pi*0.9;
R = [10; 20; 14];
v = [0.5; -0.7; 1]; 
a = [60; 20; -30]/180*pi;
trueVal = sortrows(cat(2, R, v, rad2deg(a)))

% generating code
code = hadamard(Nc);

% generating signal z
for lt = 1:size(z, 1)
    for lr = 1:size(z, 2)
        for n = 1:size(z, 3)
            for p = 1:size(z, 4)
                for t = 1:size(R, 1)
                    f_d = -2*v(t)/lambda;
                    z(lt, lr, n, p) = z(lt, lr, n, p) + exp(1j*2*pi*((2*K*R(t)/c + f_d)*mod(n - 1, 32)/f_s + 2*f_c*R(t)/c + f_d*(p - 1)*T_0 + f_c*((lt - 1)*dt + (lr - 1)*dr)*sin(a(t))/c + f_d*T_0*(lt - 1)/size(z, 1)))*code(lt, floor((n - 1)/32) + 1);
                end
            end
        end
    end
end

% receiving signal zr
zr = zeros(size(z, 2), size(z, 3), size(z, 4));
for lr = 1:size(z, 2)
    for n = 1:size(z, 3)
        for p = 1:size(z, 4)
            for lt = 1:size(z, 1)
                zr(lr, n, p) = zr(lr, n, p) + z(lt, lr, n, p);
            end
        end
    end
end
zr = awgn(zr, 10);

% distance, velocity estimation
zrr = squeeze(zr(1, :, :));
w = window2(size(zrr, 1), size(zrr, 2), @chebwin);
doppler_FFT = fft2(zrr.*w, cnt, cnt);
doppler_FFT = circshift(doppler_FFT, [0 cnt/2]).';
x = 0:rmax/cnt:rmax*(1 - 1/cnt);
y = vmax:-vmax*2/cnt:-vmax*(1 - 2/cnt);

adoppler_FFT = abs(doppler_FFT);
f1 = figure(1);
mesh(x, y, log(adoppler_FFT))
set(f1, 'Position', [10 50 500 500]);
[pks,locs] = findpeaks(log(adoppler_FFT(:)), 'SortStr', 'descend', 'NPeaks', size(R, 1), 'MinPeakHeight', 1, 'MinPeakDistance', 10^5);
rloc = zeros(size(locs)); vloc = zeros(size(locs));
rval = zeros(size(locs)); vval = zeros(size(locs));

for i = 1:size(locs)
    rloc(i) = floor(locs(i)/cnt) + 1;
    vloc(i) = mod(locs(i), cnt);
    rval(i) = x(rloc(i));
    vval(i) = y(vloc(i));
end

% angle estimation

rloc2 = zeros(size(R, 1), 1); rval2 = zeros(size(R, 1), 1);
vloc2 = zeros(size(R, 1), 1); vval2 = zeros(size(R, 1), 1);
aloc = zeros(size(R, 1), 1); aval = zeros(size(R, 1), 1);

zrr2 = zeros(size(zr, 1), Nc, size(zr, 2)/Nc, size(zr, 3));
for chp = 1:Nc
    zrr2(:, chp, :, :) = zr(:, 1 + size(zr, 2)/Nc*(chp - 1):size(zr, 2)/Nc*chp, :);
end

tx = zeros(size(z, 1), size(z, 2), size(z, 3)/Nc, size(z, 3));
tx(1, :, :, :) = (zrr2(:, 1, :, :) + zrr2(:, 2, :, :) + zrr2(:, 3, :, :) + zrr2(:, 4, :, :))/4;
tx(2, :, :, :) = (zrr2(:, 1, :, :) + zrr2(:, 2, :, :) - zrr2(:, 3, :, :) - zrr2(:, 4, :, :))/4;
tx(3, :, :, :) = (zrr2(:, 1, :, :) - zrr2(:, 2, :, :) - zrr2(:, 3, :, :) + zrr2(:, 4, :, :))/4;
tx(4, :, :, :) = (zrr2(:, 1, :, :) - zrr2(:, 2, :, :) + zrr2(:, 3, :, :) - zrr2(:, 4, :, :))/4;

zg = zeros(size(tx, 1)*size(tx, 2), size(tx, 3), size(tx, 4));
for lt = 1:size(tx, 1)
    for lr = 1:size(tx, 2)
        zg((lt - 1)*size(tx, 2) + lr, :, :) = tx(lt, lr, :, :);
    end
end
pv = zeros(size(z, 1)*size(z, 2), size(rval, 1));

% l = 1. Find locs.

zg2 = squeeze(zg(1, :, :));
w2 = window2(size(zg2, 1), size(zg2, 2), @chebwin);
d_FFT = fft2(zg2.*w2, cnt, cnt);
d_FFT = circshift(d_FFT, [0 cnt/2]).';
ad_FFT = abs(d_FFT);
[rvpks, rvlocs] = findpeaks(real(ad_FFT(:)), 'SortStr', 'descend', 'NPeaks', size(R, 1), 'MinPeakHeight', 5, 'MinPeakDistance', 10^5);
for i = 1:size(rvlocs)
    pv(1, i) = d_FFT(rvlocs(i));
    rloc2(i) = floor(rvlocs(i)/cnt) + 1;
    vloc2(i) = mod(rvlocs(i), cnt);
    rval2(i) = x(rloc2(i));
    vval2(i) = y(vloc2(i));
end

% l > 1

for l = 2:size(zg, 1)
    zg2 = squeeze(zg(l, :, :));
    d_FFT = circshift(fft2(zg2.*w2, cnt, cnt), [0 cnt/2]).';
    ad_FFT = abs(d_FFT);
    % motion compensation
    t = floor(l/size(z, 2));
    for i = 1:size(rvlocs)
        pv(l, i) = d_FFT(rvlocs(i))*exp(1j*2*pi*2*vval2(i)/lambda/size(z, 1)*T_0*t);
    end
end

f2 = figure(2);
y2 = rad2deg(asin(-1:2/(cnt - 1):1));
for i = 1:size(pv, 2)
    n_FFT = fft(pv(:, i).*chebwin(size(pv, 1)), cnt)/cnt;
    n_FFT = circshift(n_FFT, [size(n_FFT, 1)/2 0]);
    an_FFT = abs(n_FFT);
    subplot(1, size(R, 1), i);
    plot(y2, log(an_FFT));
    [apks, aloc(i)] = findpeaks(an_FFT(:), 'SortStr', 'descend', 'NPeaks', 1);
    aval(i) = y2(aloc(i));
end
set(f2, 'Position', [530 50 1000 500])

% R, V, A of targets
rv = zeros(size(rval, 1), 2);
ra = zeros(size(rval2, 1), 2);
rv(:, 1) = rval(:); rv(:, 2) = vval(:);
ra(:, 1) = rval2(:); ra(:, 2) = aval(:);
rv = sortrows(rv); ra = sortrows(ra);
measVal = cat(2, rv, ra(:, 2))

% f3 = figure(3)
% plot(x, log(abs(d_FFT)));

function w=window2(N,M,w_func)

wc=window(w_func,N);
wr=window(w_func,M);
[maskr,maskc]=meshgrid(wr,wc);
w=maskr.*maskc;

end