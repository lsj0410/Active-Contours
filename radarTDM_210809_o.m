% Motion compensated, random input. Working code

% Assume target is at 10m, 15m, 40m and moving at 1.3m/s, 0.5m/s, -0.8m/s with
% angles 5deg, -2deg, 17deg.
% 2 TX antennas, 4 RX antennas
clc; clear all;
z = zeros(2, 4, 128, 128);
c = 3*10^8;
f_c = 3*10^9; lambda = c/f_c;
B = 1.5*10^13; K = 3*10^15; f_s = 10^9;
T_c = B/K; T_0 = T_c*size(z, 1);
rmax = T_c*c/2/B*f_s; vmax = lambda/4/T_0;
dr = lambda/2; dt = dr*size(z, 2);
cnt = 1024;

% generating random input
R = rand(3, 1)*rmax;
v = (rand(3, 1) - 0.5)*vmax;
a = (rand(3, 1) - 0.5)*pi;
trueVal = sortrows(cat(2, R, v, rad2deg(a)))

for lt = 1:size(z, 1)
    for lr = 1:size(z, 2)   
        for n = 1:size(z, 3)
            for p = 1:size(z, 4)
                for i = 1:size(R, 1)   
                    f_d = -2*v(i)/lambda;
                    z(lt, lr, n, p) = z(lt, lr, n, p) + exp(1j*2*pi*((2*K*R(i)/c + f_d)*(n - 1)/f_s + 2*f_c*R(i)/c + f_d*(p - 1)*T_0 + f_c*((lt - 1)*dt + (lr - 1)*dr)*sin(a(i))/c + f_d*T_0*(lt - 1)/size(z, 1)));
                end
            end
        end
    end
end

zg = zeros(size(z, 1)*size(z, 2), size(z, 3), size(z, 4));
for lt = 1:size(z, 1)
    for lr = 1:size(z, 2)
        zg((lt - 1)*size(z, 2) + lr, :, :) = z(lt, lr, :, :);
    end
end

zg = awgn(zg, 10);

% distance, velocity estimation
zg2 = squeeze(zg(1, :, :));
w = window2(size(zg2, 1), size(zg2, 2), @chebwin);
doppler_FFT = fft2(zg2.*w);
doppler_FFT = circshift(doppler_FFT, [0 size(zg2, 2)/2]).';
x = 0:rmax/size(zg, 2):rmax*(1 - 1/size(zg, 2));
y = vmax:-vmax*2/size(zg, 3):-vmax*(1 - 2/size(zg, 3));

adoppler_FFT = abs(doppler_FFT);
f1 = figure(1);
mesh(x, y, log(adoppler_FFT))
set(f1, 'Position', [10 50 500 500]);
[pks,locs] = findpeaks(real(adoppler_FFT(:)), 'SortStr', 'descend', 'NPeaks', size(R, 1));
rloc = zeros(size(locs)); vloc = zeros(size(locs));
rval = zeros(size(locs)); vval = zeros(size(locs));

for i = 1:size(locs)
    rloc(i) = floor(locs(i)/size(zg2, 1)) + 1;
    vloc(i) = mod(locs(i), size(zg2, 1));
    rval(i) = x(rloc(i));
    vval(i) = y(vloc(i));
end

% angle estimation

tx = zeros(size(z, 1)*size(z, 2), size(rval, 1));
rloc2 = zeros(size(R, 1), 1); aloc = zeros(size(R, 1), 1);
rval2 = zeros(size(R, 1), 1); aval = zeros(size(R, 1), 1);

% l = 1. Find locs.
zg2 = squeeze(zg(1, :, :));
doppler_FFT = circshift(fft2(zg2.*w), [0 size(zg2, 2)/2]).';
adoppler_FFT = abs(doppler_FFT);
[rvpks, rvlocs] = findpeaks(real(adoppler_FFT(:)), 'SortStr', 'descend', 'NPeaks', size(R, 1));
for i = 1:size(rvlocs)
    tx(1, i) = doppler_FFT(rvlocs(i));
    rloc2(i) = floor(rvlocs(i)/size(zg2, 1)) + 1;
    rval2(i) = x(rloc2(i));
end

% l > 1
for l = 2:size(zg, 1)
    zg2 = squeeze(zg(l, :, :));
    doppler_FFT = circshift(fft2(zg2.*w), [0 size(zg2, 2)/2]).';
    adoppler_FFT = abs(doppler_FFT);
    % motion compensation
    t = floor(l/size(z, 2));
    for i = 1:size(rvlocs)
        tx(l, i) = doppler_FFT(rvlocs(i))*exp(1j*2*pi*2*vval(i)/lambda/size(z, 1)*T_0*t);
    end
end

f2 = figure(2);
y2 = rad2deg(asin(-1:2/(cnt - 1):1));
for i = 1:size(tx, 2)
    n_FFT = fft(tx(:, i), cnt)/cnt;
    n_FFT = circshift(n_FFT, [size(n_FFT, 1)/2 0]);
    an_FFT = abs(n_FFT);
    subplot(1, size(R, 1), i);
    plot(y2, an_FFT);
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

function w=window2(N,M,w_func)

wc=window(w_func,N);
wr=window(w_func,M);
[maskr,maskc]=meshgrid(wr,wc);

%maskc=repmat(wc,1,M); Old version
%maskr=repmat(wr',N,1);

w=maskr.*maskc;

end