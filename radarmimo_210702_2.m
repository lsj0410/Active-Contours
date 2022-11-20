% Attempt to separate angle graphs by object. Failed.
% Assume target is at 10m, 15m, 40m and moving at 3m/s, 1m/s, -2m/s with
% angles 5deg, -30deg, 17deg.
% 5 TX antennas, 10 RX antennas
clc; clear all;
R = [10; 15; 40];
v = [3; 1; -2.7];
a = [5 -30 17]*pi/180;
c = 3*10^8;
f_c = 3*10^9; lambda = c/f_c;
K = 3*10^15; T_0 = 5*10^-3; f_s = 10^9;
rmax = 50; vmax = 5;
z = zeros(2, 4, 100, 100);
dr = lambda/2; dt = dr*size(z, 2);

for lt = 1:size(z, 1)
    for lr = 1:size(z, 2)
        for n = 1:size(z, 3)
            for p = 1:size(z, 4)
                for i = 1:size(R)
                    f_d = -2*v(i)/lambda;
                    z(lt, lr, n, p) = z(lt, lr, n, p) + exp(1j*2*pi*(2*K*R(i)/c*(n - 1)/f_s + ((lt - 1)*dt + (lr - 1)*dr)*sin(a(i))/lambda + 2*f_c*R(i)/c + f_d*(p - 1)*T_0));
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
range_FFT = zeros(size(zg2));
doppler_FFT = zeros(size(zg2));
for i = 1:size(zg2, 1)
    range_FFT(i, :) = fft(zg2(i, :));
end
for j = 1:size(zg2, 2)
    doppler_FFT(:, j) = fft(range_FFT(:, j));
end
doppler_FFT = circshift(doppler_FFT, [0 size(zg2, 2)/2]).';
x = 0:rmax/size(zg, 2):rmax*(1 - 1/size(zg, 2));
y = vmax:-vmax*2/size(zg, 3):-vmax*(1 - 2/size(zg, 3));

adoppler_FFT = abs(doppler_FFT);
figure
mesh(x, y, adoppler_FFT)
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
        range_FFT(i, :) = fft(zg2(i, :), 100);
    end
    for j = 1:size(zg2, 2)
        doppler_FFT(:, j) = fft(range_FFT(:, j));
    end
    doppler_FFT = circshift(doppler_FFT, [0 size(zg2, 2)/2]).';
    x = 0:rmax/size(zg, 2):rmax*(1 - 1/size(zg, 2));
    y = vmax:-vmax*2/size(zg, 3):-vmax*(1 - 2/size(zg, 3));
    adoppler_FFT = abs(doppler_FFT);
    
    [pks,locs] = findpeaks(real(adoppler_FFT(:)), 'MinPeakHeight',1000, 'MinPeakProminence',1);
    rloc = zeros(size(locs)); vloc = zeros(size(locs));
    for i = 1:size(locs)
        rloc(i) = floor(locs(i)/size(zg2, 1)) + 1;
        vloc(i) = mod(locs(i), size(zg2, 1));
        zg2_2 = circshift(zg2, [0 size(zg2, 2)/2]).';
        zg1(l, i) = zg2_2(rloc(i), vloc(i));
    end
end

y2 = asin(-1:2/(size(zg1, 1) - 1):1)*180/pi;
for i = 1:size(zg1, 2)
    n_FFT = fft(zg1(:, i));
    n_FFT = circshift(n_FFT, [size(zg1, 1)/2 0]);
    an_FFT = abs(n_FFT);
    figure
    plot(y2, an_FFT);
end