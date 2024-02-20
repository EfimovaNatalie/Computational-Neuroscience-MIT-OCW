close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load fish.mat rho
%disp(rho);
num_spikes = sum(rho); % number of spikes in the trial
disp(num_spikes);
load fish.mat time
duration = time(end)*0.001; % duration of the experiment in seconds
disp(duration);
firing_rate = num_spikes/duration; % firing rate in Hz averaged over the whole experiment 
disp(firing_rate);
half_spikes = sum(rho(1:length(rho)/2)); % How many spikes are in the first half of the experiment?
%disp(half_spikes);
half_firing_rate = 2*half_spikes/duration; % What is the firing rate in Hz, averaged over the first half of the experiment? 
disp(half_firing_rate);
load fish.mat stim
max_stim = max(stim); % What is the maximum value of the stimulus
%disp(max_stim);
min_stim = min(stim); % What is the minimum value of the stimulus?
%disp(min_stim);
index = find(rho, 100); % At what time (in milliseconds) did the hundredth spike occur?
%disp(time(index(end)));
mean = sum(rho)/length(rho); % What is the mean of the spike train?
%disp(mean);
variance = var(rho); % What is the variance of the spike train? 
%disp(variance);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3rd problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len = 1000;
y = rho(1:len);
ind = time(find(rho, 160)).';
figure
plot(time(1:1000), stim(1:len), 'b')
line([ind;ind], [200; 300].*ones(size(ind)), 'Color','r')
xlabel({'time (ms)'})
ylabel({'voltage modulation'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5th problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% firing rate estimation using convolution with a boxcar filter
filter = ones(101,1)/(101); % boxcar filter
len = 10000;
signal = rho(1:len);
est_firing_rate = conv(filter, signal);
prob = est_firing_rate(51:end-50);
dt = time(end)/length(time);
frequency_coeff = 1/(dt*0.001); % 0.001 to get Hz
rate = prob*frequency_coeff;
%disp(length(prob));
figure
plot(time(1:len), rate, 'b', time(1:len), stim(1:len), 'r')
legend({'firing rate','stimulus'})
xlabel('time')
ylabel({'voltage modulation / firing rate'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6th problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = polyfit(stim(1:len), rate, 1); % Using the polyfit function, find the coefficients a and b
rate_approx = polyval(p, stim(1:len)); % a*stim+b best approximates prob
figure
plot(time(1:len), rate, 'r', time(1:len), rate_approx, 'b') % a*stim + b and est firing rate comparison
legend({'stimulus', 'firing rate', 'approx firing rate'})
xlabel({'time (ms)'})
ylabel({'firing rate (Hz)'})
figure
plot(stim(1:len), rate_approx, 'LineWidth', 2) % straight line fitted to the data points.
hold on
s = scatter(stim(1:len), rate, "filled"); % data points 
s.SizeData = 5;
s.MarkerFaceColor = 'g';
legend({'a*stim + b', 'data points'})
xlabel({'voltage modulation'})
ylabel({'firing rate (Hz)'})
R = corrcoef(stim(1:len), rate); % correlation coefficient calculating
title('correlation coeff = ', R(2));
disp(R(2));
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparison of different window lengths for filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filter_2 = ones(1001,1)/(1001); % boxcar filter
est_firing_rate_2 = conv(filter_2, signal);
prob_2 = est_firing_rate_2(501:end-500);
rate_2 = prob_2*frequency_coeff;
R_2 = corrcoef(stim(1:len), rate_2); % correlation coefficient calculating
disp(R_2(2)); 

filter_3 = ones(11,1)/(11); % boxcar filter
est_firing_rate_3 = conv(filter_3, signal);
prob_3 = est_firing_rate_3(6:end-5);
rate_3 = prob_3*frequency_coeff;
R_3 = corrcoef(stim(1:len), rate_3); % correlation coefficient calculating
disp(R_3(2)); 

figure
plot(time(1:len), rate, 'b', time(1:len), rate_2, time(1:len), rate_3, time(1:len), stim(1:len), 'r')
legend({'101', '1001', '11','stimulus'})
xlabel('time')
ylabel({'voltage modulation / firing rate'})
% We can find that filter with window length = 1001 gives the best fit (corr coeff is maximum, =~ 0.83)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7th problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_rho = polyfit(stim(1:len), rho(1:len), 1);
rho_approx = polyval(p_rho, stim(1:len)); % p contains a and b, a*stim+b best approximates rho
figure
plot(stim(1:len), rho_approx, 'r', 'LineWidth', 2)
hold on
s = scatter(stim(1:len), rho(1:len), 'filled'); 
s.SizeData = 5;
s.MarkerFaceColor = 'b';
legend({'a*stim + b', 'data points'})
xlabel('stimulus')
ylabel('rho')
hold off
%R_rho = corrcoef(stim(1:len), rho(1:len)); % fit is bad
%disp(R_rho(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8th problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 1000; % max lag to compute autocorrelation
autocorr_stim = xcorr(stim, stim, M); % autocorrelation of the stimulus
figure
plot(autocorr_stim) % central peak corresponds to zero lag, its width ~ 300 ms
title('Stimulus Autocorrelation')
autocorr_rho = xcorr(rho, rho, M);
figure
plot(autocorr_rho) % we can find much more peaks than for stimulus due to binary data and quite monotonous (similar peak's height and often appear) signal
title('Spike Train Aturocorrelation') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9th problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cross_corr = xcorr(stim, rho, M); % cross-correlation of the stimulus and spike train
figure
plot(cross_corr)
title('Cross-correlation') 