function psd_resynthesis(coeffs_in,fs,cin)
%coeffs_in(coeffs_in<0) = 0;
coeffs_linear = 10.^(coeffs_in);
samp_rate = 10000;
sig_length = 1;
t =0:1./samp_rate:(sig_length-1./samp_rate);

sum_signal = zeros(1,numel(t));
sum_signal_log = zeros(1,numel(t));

for mm = fs
sum_signal = sum_signal+sqrt(coeffs_linear(mm))*abs(sin(mm*t*2*pi ));
sum_signal_log = sum_signal_log+(coeffs_in(mm))./2*(sin(mm*t*2*pi ));

end

figure(99)
hold on
plot(t,sum_signal,cin)

figure(100)
hold on
plot(t,sum_signal_log,cin)