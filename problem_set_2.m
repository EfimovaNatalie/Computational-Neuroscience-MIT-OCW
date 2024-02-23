close all
load fish.mat rho
load fish.mat stim
load fish.mat time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% task 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 10; % vectors' length
g = rand(m, 1);
h = rand(m, 1);

corr = my_xcorr(g, h);

A = cat(1, corr', xcorr(g, h)'); % matrix for comparison xcorr and my_xcorr functions' performance

%disp(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% task 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_lag = 300;
mean_stim = sum(stim)/length(stim);
y = stim(1:1000) - mean_stim*ones(1000, 1);
mean_sp_train = sum(rho)/length(rho);
x = rho(1:1000) - mean_sp_train*ones(1000, 1);

cross_corr = xcorr(x, y, max_lag);
filter_h = cross_corr(200:end)/sum(rho(1:1000));
figure
plot(filter_h); % filter plot corresponding to the indices range from M1 = -100 to M2 = 300
ylabel('filter h')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% task 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
approx_y = conv(x, filter_h);
%disp(length(approx_y))
figure
plot(time(1:1000), approx_y(67:end-335), time(1:1000), y) % (h * x) convolution approximates y
xlabel('time (ms)')
ylabel('voltage modulation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% task 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_xx = xcorr(x, x); % computing auto-correlation
C_xy = xcorr(x, y); % computing cross-correlation
C_xx = toeplitz(c_xx);
wiener_filter_h = (C_xy)\(C_xx); % defining filter h from Wiener-Hopf equations
figure
plot(wiener_filter_h)
ylabel('Wiener filter h')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% task 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
approx_w_y = conv(x, wiener_filter_h, 'same');
figure
plot(time(1:1000), approx_w_y, time(1:1000), y)
xlabel('time (ms)')
ylabel('voltage modulation')

function correlation = my_xcorr(g, h) % my own xcorr function using conv command
    correlation = flip(conv(h, flip(g)));
end