
function n = findpeaks(x)
% Find peaks.
% n = findpeaks(x)
diff_x = diff(x);

diff_x_greater_than_zero = diff_x > 0;          % increasing, positive gradient

diff_diff = diff(diff_x_greater_than_zero);         % find on and off points for the positive gradient

n    = find(diff_diff < 0);                               % this is the samples with the gradient goes from +ve to negative. maxima

u    = find(x(n+1) > x(n));                     % find where x(n) is increasing, I think we could replace these two lines with n=n+1;

n(u) = n(u)+1;
% 
% plot(diff_x_greater_than_zero)
% hold
% plot(diff_diff,'r')
% xlabel('Samples')
% legend('Samples with +ve gradient','On off times of +ve gradient')
% ylim([-1.1 1.1])
