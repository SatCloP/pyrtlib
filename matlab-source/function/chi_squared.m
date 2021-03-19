function [x2,tb1,tb2,expected,observed] = chi_squared(samples);

% function [x2,tb1,tb2,expected,observed] = chi_squared(samples);
%
% returns the chi-squared value x2 for a Gaussian distribution of samples,
% along with the expected and observed number of occurances for the std 
% bins defined by tb1 and tb2.

% Probability p that a measurement fall within t standard deviations of 
% the mean for a true Gaussian distribution
t = [0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00 2.50 3.00 3.50  4.00]';
p = [0.00 20   38   55   68   79   87   92   96.4 98.8 99.7 99.95 99.99]';

plot(t,p,'o')

t = [-flipud(t(2:length(t)));t];
tb1 = [-inf;t]; % create bin boundaries
tb2 =  [t;inf];
pbin=[diff(p)]; % probability that a measurement fall only w/in a bin
pbin = [flipud(pbin);pbin];
pbin = [100-99.99; pbin; 100-99.9];

N = length(samples);
mn = mean(samples);
st = std(samples);

nbins = length(pbin);

disp(['number of bins = ' num2str(nbins)])

expected = pbin/100*N;   % expected number of occurances

observed = zeros(nbins,1);
for i = 1:nbins
observed(i) = length(find(samples > mn+tb1(i)*st & samples <= mn+tb2(i)*st));
end

% plot(t,expected(1:nbins-1),t,observed(1:nbins-1));

x2 = sum(((observed-expected).^2)./expected);

% keyboard

