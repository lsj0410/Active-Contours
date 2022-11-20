% Working code. Outputs d-v graphs, angle graphs of each object and angle
% values. Not logscale. Does not pair peaks.
% Assume target is at 10m, 15m, 40m and moving at 3m/s, 1m/s, -2m/s with
% angles 5deg, -30deg, 17deg.
% 2 TX antennas, 8 RX antennas
clc; clear all;
zg = zeros(8, 100, 100);
R = [10; 15; 40];
v = [2.5; 1; -1.2];
a = [50 -30 17]*pi/180;
c = 3*10^8;
f_c = 3*10^9; lambda = c/f_c;
B = 1.5*10^13; K = 3*10^15; f_s = 10^9;
T_c = B/K;
rmax = T_c*c/2/B*f_s; vmax = lambda/4/T_c;
dr = lambda/2;

for l = 1:size(zg, 1)
    for n = 1:size(zg, 2)
        for p = 1:size(zg, 3)
            for i = 1:size(R, 1)
                f_d = -2*v(i)/lambda;
                zg(l, n, p) = zg(l, n, p) + exp(1j*2*pi*((2*K*R(i)/c + f_d)*(n - 1)/f_s + 2*f_c*R(i)/c + f_d*(p - 1)*T_c + f_c*(l - 1)*dr*sin(a(i))/c));
            end
        end
    end
end

zg = awgn(zg, 10);

% distance, velocity estimation
zg2 = squeeze(zg(1, :, :));
w2 = window2(size(zg2, 1), size(zg2, 2), @chebwin);
doppler_FFT = circshift(fft2(zg2.*w2), [0 size(zg2, 1)/2]).';
x = 0:rmax/size(zg, 2):rmax*(1 - 1/size(zg, 2));
y = vmax:-vmax*2/size(zg, 3):-vmax*(1 - 2/size(zg, 3));

adoppler_FFT = abs(doppler_FFT);
f1 = figure(1);
mesh(x, y, log(adoppler_FFT))
set(f1, 'Position', [10 10 500 500])
[pks,locs] = findpeaks(real(adoppler_FFT(:)), 'MinPeakHeight',1000, 'MinPeakProminence',1);
rloc = zeros(size(locs)); vloc = zeros(size(locs));
rval = zeros(size(locs)); vval = zeros(size(locs));

for i = 1:size(locs)
    rloc(i) = floor(locs(i)/size(zg2, 1)) + 1;
    vloc(i) = mod(locs(i), size(zg2, 1));
    rval(i) = x(rloc(i));
    vval(i) = y(vloc(i));
end

% angle estimation
zg1 = zeros(size(zg, 1), size(rval, 1));
for l = 1:size(zg, 1)
    zg2 = squeeze(zg(l, :, :));
    range_FFT = zeros(size(zg2));
    doppler_FFT = zeros(size(zg2));
    for i = 1:size(zg2, 1)
        range_FFT(i, :) = fft(zg2(i, :));
    end
    for j = 1:size(zg2, 2)
        doppler_FFT(:, j) = fft(range_FFT(:, j));
    end
    adoppler_FFT = abs(doppler_FFT);
    [pks,locs] = findpeaks(real(adoppler_FFT(:)), 'MinPeakHeight',1000, 'MinPeakProminence',1);
    
    for i = 1:size(locs)
        zg1(l, i) = doppler_FFT(locs(i));
    end
end

y2 = asin(-1:2/(100 - 1):1)*180/pi;
f2 = figure(2);
for i = 1:size(zg1, 2)
    n_FFT = fft(zg1(:, i).*chebwin(8), 100);
    n_FFT = circshift(n_FFT, [50 0]);
    an_FFT = abs(n_FFT);
    subplot(1, size(R, 1), i);
    plot(y2, log(an_FFT));
    [pks,locs] = findpeaks(real(an_FFT(:)), 'SortStr', 'descend', 'NPeaks', 1);
    y2(locs)
end
set(f2, 'Position', [550 50 1000 1000/size(R, 1)]);

function w=window2(N,M,w_func)

wc=window(w_func,N);
wr=window(w_func,M);
[maskr,maskc]=meshgrid(wr,wc);
w=maskr.*maskc;

end