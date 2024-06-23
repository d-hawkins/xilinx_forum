% -----------------------------------------------------------------------------
% quad_core_adc_spectra.m
%
% 6/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% Quad core ADC spectra.
%
% This script analyzes ideal data for an ADC constructed by interlacing
% four ADC cores. Each core has matched gain, and the core clocks are
% skewed by exactly one quarter of the per core clock period.
%
% The raw ADC samples are creating by filtering real-valued Gaussian noise,
% at the interlaced sampling rate, and then decimate-by-4 to generate the
% per ADC samples.
%
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
% Example
% -----------------------------------------------------------------------------
%
% Example number per the github README.md description.
%
% Example   Description
% -------   -----------
%    1      Four cores measuring 4x the bandwidth of a single core
%    2      Four cores oversampling the first Nyquist zone
%    3      Four cores oversampling the second Nyquist zone
%    4      Four cores oversampling the third Nyquist zone
%    5      Four cores oversampling the forth Nyquist zone
%    6      Four cores oversampling the first and second Nyquist zone
%    7      Four cores oversampling the first and third Nyquist zone
%
example = 7;

% Generate figures for the README.md file
generate_figures = 1;

% -----------------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------------
%
% Number of ADC cores
Nc = 4;

% Number of time samples per spectra (per core)
Nt = 1024;

% Number of estimates to incoherently average (for spectrum calculation)
Ne = 1000;

% Bit-width
Bx = 16;

% -----------------------------------------------------------------------------
% Power spectra
% -----------------------------------------------------------------------------
%
% Input noise loading factor (LF)
LF_dB = -15;

% Loading factor (normalized RMS)
LF = 10^(LF_dB/20);

% Uniform (white) or bandpass spectra?
use_bandpass = 1;

% -----------------------------------------------------------------------------
% FFT window function
% -----------------------------------------------------------------------------
%
% The frequency channel response of the FFT is sinc-like. This sinc-like
% response has such high sidelobes, that the stopband of the filtered input
% noise signal gets filled in, so you cannot see the quantization noise floor.
% This m-file windows the floating-point data and then quantizes it. This
% reduces the sidelobe response, but leaves the quantization noise floor where
% it is expected.
%

% Window function
w = kaiser(Nt+1,25);
w = w(1:Nt)';      % DFT cyclic symmetry

% Remove the conditional condition to see the unwindowed sinc
% response spectrum with the bandpass spectra
if (use_bandpass ~= 1)
	w = ones(1,Nt);
end

% Power spectrum normalization
w_cg    = mean(w);          % Coherent gain
w_ig    = sqrt(mean(w.^2)); % Incoherent gain
w_cg_dB = 20*log10(w_cg);
w_ig_dB = 20*log10(w_ig);

% -----------------------------------------------------------------------------
% Filter
% -----------------------------------------------------------------------------
%
% Create a filter passband within each ADC cores Nyquist zone.
%
% Number of coefficients in the filter
Nh = 4096;

% Filter coefficients response window
wh = kaiser(Nh,30).';

% Bandwidth per core ADC Nyquist zone
Bh = 2*round(0.8*Nh/(4*Nc));

% Slope over the band
%SdB = [0  0];
SdB = [0 20];

% -250dB out-of-band response
H_dB = -250*ones(1,Nh);

% Positive frequencies
% * Each core's Nyquist zone is Nh/(2*Nc) samples
% * Center of the first band is Nh/(4*Nc)
n_offset = Nh/(4*Nc);
n_step = Nh/(2*Nc);
n_span = [-Bh/2:Bh/2];
for n = 1:Nc,

	% Skip bands
	if (example == 2)
		if ((n==2) || (n==3) || (n==4))
			continue
		end
	elseif (example == 3)
		if ((n==1) || (n==3) || (n==4))
			continue
		end
	elseif (example == 4)
		if ((n==1) || (n==2) || (n==4))
			continue
		end
	elseif (example == 5)
		if ((n==1) || (n==2) || (n==3))
			continue
		end
	elseif (example == 6)
		if ((n==3) || (n==4))
			continue
		end
	elseif (example == 7)
		if ((n==2) || (n==4))
			continue
		end
	end

	H_dB(1+n_offset+(n-1)*n_step+n_span) = SdB(1) + [0:Bh]/Bh*(SdB(2)-SdB(1));
end

% Copy positive to negative frequencies
H_dB((Nh/2+2):Nh) = H_dB(Nh/2:-1:2);

% Regular frequency format
h_channel = [-Nh/2:Nh/2-1];
H_dB = fftshift(H_dB);

% Plot the filter ideal response
figure(1)
clf
for n = 1:2*Nc-1,
	plot(h_channel(n*n_step+1)*[1 1], [-300 50], 'k--')
	if (n == 1)
		hold on
	end
end
plot(h_channel,H_dB,'b')
xlabel('Frequency Channel')
ylabel('Power (dB)')
axis([-inf inf -inf inf])
title('Filter Response')

% Linear response
H = 10.^(H_dB/20);

% Windowed coefficients response
h = wh.*fftshift(ifft(fftshift(H)));

% Power spectra
H = fftshift(fft(fftshift(h)));
h_ig = sqrt(mean(abs(H).^2));
H_dB = 10*log10(abs(H).^2);

% Plot the filter actual response
figure(1)
plot(h_channel,H_dB,'r')

% -----------------------------------------------------------------------------
% Filtered Signal
% -----------------------------------------------------------------------------
%
fprintf('Create the input signal (filtered noise)\n');

% Random number generator
rng(1234,"twister");

% Total number of samples to filter
N = Nt*Ne + Nh;

% Real-valued noise
x = randn(1,N);

% Filter
% * comment out this filtering step to see the effect
%   of quantization on wideband noise
if (use_bandpass)
	x = filter(h,1,x);
end

% Discard the transient
x = x(Nh + [1:Nt*Ne]);

% Scale
x_std = std(x);
x     = LF*(2^(Bx-1)-1)*x/x_std;

% Apply the window function to FFT data blocks BEFORE quantization
%
% The window function has to be applied before the data is quantized,
% so that if data is passed to a fixed-point FFT, quantization does not
% occur *twice* (which raises the quantization noise floor by 6dB).
%
% The application of this window reduces the signal power by the incoherent
% gain of the window. You cannot compensate the input signal for the loss
% of this power, since increasing the amplitude causes clipping.
%
for n = 1:Ne,
	x( (n-1)*Nt + [1:Nt] ) = w.*x( (n-1)*Nt + [1:Nt] );
end

% Save a copy of the unquantized but windowed data
xf = x;

% Quantize/saturate
x = round(x);
m = find(x > 2^(Bx-1)-1);
if ~isempty(m)
	fprintf('WARNING: saturated %d samples to maximum\n', length(m));
	x(m) = 2^(Bx-1)-1;
end
m = find(x < -2^(Bx-1)+1);
if ~isempty(m)
	fprintf('WARNING: saturated %d samples to minimum\n', length(m));
	x(m) = -2^(Bx-1)+1;
end

% Quantization noise
%  x = xf + xq
xq = x - xf;

% =============================================================================
% Quad-Core ADC Spectra
% =============================================================================
%
% Calculate the spectra for all cores interlaced.
%
% -----------------------------------------------------------------------------
% Power Spectrum Calculation
% -----------------------------------------------------------------------------
%
fprintf('Calculate the power spectra\n');

Rxx = zeros(1,Nt);
Rff = zeros(1,Nt);
Rqq = zeros(1,Nt);

for n = 1:Ne,
	% MATLAB FFT
	x_est = x( (n-1)*Nt + [1:Nt] );
	X_est = fftshift(fft(x_est));
	Rxx = Rxx + abs(X_est).^2;

	% Full-precision (no quantization)
	f_est = xf( (n-1)*Nt + [1:Nt] );
	F_est = fftshift(fft(f_est));
	Rff = Rff + abs(F_est).^2;

	% Quantization noise
	q_est = xq( (n-1)*Nt + [1:Nt] );
	Q_est = fftshift(fft(q_est));
	Rqq = Rqq + abs(Q_est).^2;

	% Progress indication
	fprintf('.');
	if ( mod(n,50) == 0)
		fprintf(' %d\n', n);
	end
end
fprintf('\n')

% Normalize
%
% * Ne for the number of estimates
% * (2^(Bx-1))^2 to scale to fractional integer
% * Nt to normalize noise
%
Rxx = Rxx/(Ne*Nt*(2^(Bx-1))^2);
Rff = Rff/(Ne*Nt*(2^(Bx-1))^2);
Rqq = Rqq/(Ne*Nt*(2^(Bx-1))^2);

% dB with -200dB noise floor
Rxx_dB = 10*log10(Rxx + 10^(-20));
Rff_dB = 10*log10(Rff + 10^(-20));
Rqq_dB = 10*log10(Rqq + 10^(-20));

% -----------------------------------------------------------------
% Figures
% -----------------------------------------------------------------
%
n = [-Nt/2:Nt/2-1];

% Starting figure number
%  * change 'fignum' when comparing different simulation runs
fignum = 1;

% --------------------------
% Power Spectra
% --------------------------
%
figure(fignum+1)
clf
% 0dB axis
plot(Nt/2*[-1 1], [0 0], 'k--')
hold on

ph(1) = plot(n,Rff_dB,'g');
ph(2) = plot(n,Rxx_dB,'b');
ph(3) = plot(n,Rqq_dB,'r');

% Noise signal level
% (reduced by the window incoherent gain)
plot(Nt/2*[-1 1], [1 1]*(LF_dB+w_ig_dB), 'b--')

% Quantization noise floor
QN_dB = -6.02*Bx - 4.77;
if (use_bandpass == 1)
	% Adjust down by 1dB when the signal is not uniform bandpass
	QN_dB = QN_dB - 1;
end
plot(Nt/2*[-1 1], [1 1]*QN_dB, 'r--')

y_min = floor(QN_dB/10)*10;
axis([-Nt/2 Nt/2 y_min 20])

xlabel('Frequency Channel')
ylabel('Loading Factor (dB)')
title('Power spectra for all cores interlaced')

hl=legend(ph, ...
	'Floating-point','Fixed-point','Quantization', ...
	'Location','North','AutoUpdate','Off');
set(hl,'FontSize',10);

% Annotate the figure
%
% * make sure that any axis calls (which re-scale the figure)
%   occur prior to any call to ds2nfu, since the rescaling
%   changes the normalized figure units.
%
% Input signal power
str = sprintf('%.2fdB',LF_dB+w_ig_dB);
text(0.95*(-Nt/2),LF_dB+w_ig_dB+2,str, ...
	'HorizontalAlignment','left', ...
	'VerticalAlignment','bottom', ...
	'FontSize',10, ...
	'Color','b')

% Quantization noise floor
str = sprintf('%.2fdB',QN_dB);
text(0.95*(-Nt/2),QN_dB+2,str, ...
	'HorizontalAlignment','left', ...
	'VerticalAlignment','bottom', ...
	'FontSize',10, ...
	'Color','r')

if (generate_figures)
	if (use_bandpass)
		pngname = sprintf('ex%d_bandpass_spectra.png', example);
	else
		pngname = sprintf('ex%d_uniform_spectra.png', example);
	end
	exportgraphics(gcf,pngname,'Resolution',300)
end

% =============================================================================
% Single-Core ADC Spectra
% =============================================================================
%
% -----------------------------------------------------------------------------
% Decimated samples
% -----------------------------------------------------------------------------
%
yf = zeros(Nc,Ne*Nt/Nc);
y  = zeros(Nc,Ne*Nt/Nc);
yq = zeros(Nc,Ne*Nt/Nc);
for n = 1:Nc,
	yf(n,:) = xf((n-1)+[1:Nc:Ne*Nt]);
	y(n,:)  =  x((n-1)+[1:Nc:Ne*Nt]);
	yq(n,:) = xq((n-1)+[1:Nc:Ne*Nt]);
end

% -----------------------------------------------------------------------------
% Power Spectrum Calculation
% -----------------------------------------------------------------------------
%
fprintf('Calculate the power spectra\n');

Dxx = zeros(4,Nt/Nc);
Dff = zeros(4,Nt/Nc);
Dqq = zeros(4,Nt/Nc);

for m = 1:Nc,
	for n = 1:Ne,
		% MATLAB FFT
		x_est = y(m, (n-1)*Nt/Nc + [1:Nt/Nc] );
		X_est = fftshift(fft(x_est));
		Dxx(m,:) = Dxx(m,:) + abs(X_est).^2;

		% Full-precision (no quantization)
		f_est = yf(m, (n-1)*Nt/Nc + [1:Nt/Nc] );
		F_est = fftshift(fft(f_est));
		Dff(m,:) = Dff(m,:) + abs(F_est).^2;

		% Quantization noise
		q_est = yq(m, (n-1)*Nt/Nc + [1:Nt/Nc] );
		Q_est = fftshift(fft(q_est));
		Dqq(m,:) = Dqq(m,:) + abs(Q_est).^2;

		% Progress indication
		fprintf('.');
		if ( mod(n,50) == 0)
			fprintf(' %d\n', n);
		end
	end
end
fprintf('\n')

% Normalize
%
% * Ne for the number of estimates
% * (2^(Bx-1))^2 to scale to fractional integer
% * Nt/Nc to normalize noise
%
Dxx = Dxx/(Ne*Nt/Nc*(2^(Bx-1))^2);
Dff = Dff/(Ne*Nt/Nc*(2^(Bx-1))^2);
Dqq = Dqq/(Ne*Nt/Nc*(2^(Bx-1))^2);

% dB with -200dB noise floor
Dxx_dB = 10*log10(Dxx + 10^(-20));
Dff_dB = 10*log10(Dff + 10^(-20));
Dqq_dB = 10*log10(Dqq + 10^(-20));

% -----------------------------------------------------------------
% Figures
% -----------------------------------------------------------------
%
n = [-Nt/(2*Nc):Nt/(2*Nc)-1];

% Starting figure number
%  * change 'fignum' when comparing different simulation runs
fignum = 2;

% --------------------------
% Power Spectra
% --------------------------
%
for m = 1:Nc,
	fignum = fignum+1;
	figure(fignum)
	clf
	% 0dB axis
	plot(Nt/2*[-1 1], [0 0], 'k--')
	hold on

	ph(1) = plot(n,Dff_dB(m,:),'g');
	ph(2) = plot(n,Dxx_dB(m,:),'b');
	ph(3) = plot(n,Dqq_dB(m,:),'r');

	% Noise signal level
	% (reduced by the window incoherent gain)
	plot(Nt/(2*Nc)*[-1 1], [1 1]*(LF_dB+w_ig_dB), 'b--')

	% Quantization noise floor
	QN_dB = -6.02*Bx - 4.77;
	if (use_bandpass == 1)
		% Adjust down by 1dB when the signal is not uniform bandpass
		QN_dB = QN_dB - 1;
	end
	plot(Nt/(2*Nc)*[-1 1], [1 1]*QN_dB, 'r--')

	y_min = floor(QN_dB/10)*10;
	axis([-Nt/(2*Nc) Nt/(2*Nc) y_min 20])

	xlabel('Frequency Channel')
	ylabel('Loading Factor (dB)')
	title(sprintf('Power spectra for core %d', m))

	hl=legend(ph, ...
		'Floating-point','Fixed-point','Quantization', ...
		'Location','North','AutoUpdate','Off');
	set(hl,'FontSize',10);

	% Annotate the figure
	%
	% * make sure that any axis calls (which re-scale the figure)
	%   occur prior to any call to ds2nfu, since the rescaling
	%   changes the normalized figure units.
	%
	% Input signal power
	str = sprintf('%.2fdB',LF_dB+w_ig_dB);
	text(0.95*(-Nt/(2*Nc)),LF_dB+w_ig_dB+2,str, ...
		'HorizontalAlignment','left', ...
		'VerticalAlignment','bottom', ...
		'FontSize',10, ...
		'Color','b')

	% Quantization noise floor
	str = sprintf('%.2fdB',QN_dB);
	text(0.95*(-Nt/(2*Nc)),QN_dB+2,str, ...
		'HorizontalAlignment','left', ...
		'VerticalAlignment','bottom', ...
		'FontSize',10, ...
		'Color','r')

	% Just print core#1, since they all look the same
	if (m~=1)
		continue
	end
	if (generate_figures)
		if (use_bandpass)
			pngname = sprintf('ex%d_bandpass_spectra_core%d.png', example, m);
		else
			pngname = sprintf('ex%d_uniform_spectra_core%d.png', example, m);
		end
		exportgraphics(gcf,pngname,'Resolution',300)
	end
end

% =============================================================================
% IQ ADC Spectra
% =============================================================================
%
% Treat the data from the two of the ADCs as I and Q
%
% https://support.xilinx.com/s/question/0D54U00008UxFy9SAF/we-are-using-the-logicore-jesd204-core-to-interface-to-adcs-we-have-a-phase-delay-between-the-2-adcs-to-support-iq-sampling-can-the-core-manage-a-phase-offset-between-the-2-channels-of-this-much?language=en_US
%
% Using cores 1+2 or cores 1+3 does not work (as expected).
%
% Using cores 1+3 to interlace real-valued samples at twice the per core rate
% would work.
%
if (1)
	% Cores 90-degrees out-of-phase
	zf = yf(1,:) + 1j*yf(2,:);
	z  =  y(1,:) + 1j* y(2,:);
	zq = yq(1,:) + 1j*yq(2,:);
else
	% Cores 180-degrees out-of-phase
	zf = yf(1,:) + 1j*yf(3,:);
	z  =  y(1,:) + 1j* y(3,:);
	zq = yq(1,:) + 1j*yq(3,:);
end

% -----------------------------------------------------------------------------
% Power Spectrum Calculation
% -----------------------------------------------------------------------------
%
fprintf('Calculate the power spectra\n');

Zxx = zeros(4,Nt/Nc);
Zff = zeros(4,Nt/Nc);
Zqq = zeros(4,Nt/Nc);

for n = 1:Ne,
	% MATLAB FFT
	x_est = z((n-1)*Nt/Nc + [1:Nt/Nc] );
	X_est = fftshift(fft(x_est));
	Zxx = Zxx + abs(X_est).^2;

	% Full-precision (no quantization)
	f_est = zf((n-1)*Nt/Nc + [1:Nt/Nc] );
	F_est = fftshift(fft(f_est));
	Zff = Zff + abs(F_est).^2;

	% Quantization noise
	q_est = zq((n-1)*Nt/Nc + [1:Nt/Nc] );
	Q_est = fftshift(fft(q_est));
	Zqq = Zqq + abs(Q_est).^2;

	% Progress indication
	fprintf('.');
	if ( mod(n,50) == 0)
		fprintf(' %d\n', n);
	end
end
fprintf('\n')

% Normalize
%
% * Ne for the number of estimates
% * (2^(Bx-1))^2 to scale to fractional integer
% * Nt/Nc to normalize noise
%
Zxx = Zxx/(Ne*Nt/Nc*(2^(Bx-1))^2);
Zff = Zff/(Ne*Nt/Nc*(2^(Bx-1))^2);
Zqq = Zqq/(Ne*Nt/Nc*(2^(Bx-1))^2);

% dB with -200dB noise floor
Zxx_dB = 10*log10(Zxx + 10^(-20));
Zff_dB = 10*log10(Zff + 10^(-20));
Zqq_dB = 10*log10(Zqq + 10^(-20));

% -----------------------------------------------------------------
% Figures
% -----------------------------------------------------------------
%
n = [-Nt/(2*Nc):Nt/(2*Nc)-1];

% Starting figure number
%  * change 'fignum' when comparing different simulation runs
fignum = 6;

% --------------------------
% Power Spectra
% --------------------------
%
fignum = fignum+1;
figure(fignum)
clf
% 0dB axis
plot(Nt/2*[-1 1], [0 0], 'k--')
hold on

ph(1) = plot(n,Dff_dB(m,:),'g');
ph(2) = plot(n,Dxx_dB(m,:),'b');
ph(3) = plot(n,Dqq_dB(m,:),'r');

% Noise signal level
% (reduced by the window incoherent gain)
plot(Nt/(2*Nc)*[-1 1], [1 1]*(LF_dB+w_ig_dB), 'b--')

% Quantization noise floor
QN_dB = -6.02*Bx - 4.77;
if (use_bandpass == 1)
	% Adjust down by 1dB when the signal is not uniform bandpass
	QN_dB = QN_dB - 1;
end
plot(Nt/(2*Nc)*[-1 1], [1 1]*QN_dB, 'r--')

y_min = floor(QN_dB/10)*10;
axis([-Nt/(2*Nc) Nt/(2*Nc) y_min 20])

xlabel('Frequency Channel')
ylabel('Loading Factor (dB)')
title('Power spectra for two cores processed as I+Q')

hl=legend(ph, ...
	'Floating-point','Fixed-point','Quantization', ...
	'Location','North','AutoUpdate','Off');
set(hl,'FontSize',10);

% Annotate the figure
%
% * make sure that any axis calls (which re-scale the figure)
%   occur prior to any call to ds2nfu, since the rescaling
%   changes the normalized figure units.
%
% Input signal power
str = sprintf('%.2fdB',LF_dB+w_ig_dB);
text(0.95*(-Nt/(2*Nc)),LF_dB+w_ig_dB+2,str, ...
	'HorizontalAlignment','left', ...
	'VerticalAlignment','bottom', ...
	'FontSize',10, ...
	'Color','b')

% Quantization noise floor
str = sprintf('%.2fdB',QN_dB);
text(0.95*(-Nt/(2*Nc)),QN_dB+2,str, ...
	'HorizontalAlignment','left', ...
	'VerticalAlignment','bottom', ...
	'FontSize',10, ...
	'Color','r')

if (generate_figures)
	if (use_bandpass)
		pngname = sprintf('ex%d_bandpass_spectra_iq.png', example);
	else
		pngname = sprintf('ex%d_uniform_spectra_iq.png', example);
	end
	exportgraphics(gcf,pngname,'Resolution',300)
end
