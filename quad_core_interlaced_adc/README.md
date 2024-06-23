# Quad-core Interlaced ADC Example

6/22/2024 D. W. Hawkins (dwh@caltech.edu)

## Introduction

The Xilinx forum discussion:

"We are using the logiCORE JESD204 core to interface to ADCs. We have a phase delay between the 2 ADCs to support IQ sampling. Can the core manage a phase offset between the 2 channels of this much?"

https://support.xilinx.com/s/question/0D54U00008UxFy9SAF/we-are-using-the-logicore-jesd204-core-to-interface-to-adcs-we-have-a-phase-delay-between-the-2-adcs-to-support-iq-sampling-can-the-core-manage-a-phase-offset-between-the-2-channels-of-this-much

Indicated they were trying to use a 1/4 delay between two ADC clocks sampling real-valued signals to perform IQ demodulation.

If each ADC was sampling real-valued wideband noise at F_core, then the maximum bandwidth per core sampled would be F_core/2. 

If two cores were interlaced, with the second device clocked 1/2 period after the first, then the real-valued samples could be interlaced-by-two to effectively sample the input signal at 2*F_core, with a maximum sampled bandwidth of F_core.

If four cores were interlaced, with each device clocked 1/4 period after the previous device, then the real-valued samples could be interlaced-by-four to effectively sample the input signal at 4*F_core, with a maximum sampled bandwidth of 2*F_core.

The MATLAB script in this directory demonstrates this.

This README.md can be viewed in a web browser looking at the github page, or a local version can be viewed using a Markdown viewer, eg., in Visual Studio Code use "ctrl+k, v".

## Examples

### Example 1: Quad-core Spectra

The following figure shows the power spectra from Gaussian noise samples that were filtered to create four separate bands within four Nyquist zones. This is the ideal data for a quad-core interlaced ADC.

<p align="center">
  <img src="png/ex1_bandpass_spectra.png" alt="ex1_bandpass_spectra" style="width:600px;"/>
</p>

The following figure shows the power spectra calculated from the data for one core, i.e., the data generated for the quad-core ADC was decimated-by-4 for each of the cores, and this is the spectra for core 1 (the spectra for cores 2, 3, and 4) look similar. The "bathtub" spectra is created because Nyquist zones 1, 2, 3, and 4 alias on top of each other, with Nyquist zones 2 and 4 undergoing frequency reversal during the aliasing.

<p align="center">
  <img src="png/ex1_bandpass_spectra_core1.png" alt="ex1_bandpass_spectra_core1" style="width:600px;"/>
</p>

The following figure shows the power spectra for core 1 and core 2 processed as IQ data. The power spectra shows the same aliases response as the core 1 spectra.

<p align="center">
  <img src="png/ex1_bandpass_spectra_iq.png" alt="ex1_bandpass_spectra_iq" style="width:600px;"/>
</p>

---
### Example 2: Oversampled first Nyquist zone

In this example, the quad-core ADC oversamples a signal in the first Nyquist zone. The core 1 and IQ spectra look similar (but that does not mean the IQ spectra is correct!).

<p align="center">
  <img src="png/ex2_bandpass_spectra.png"       alt="ex2_bandpass_spectra" style="width:600px;"/>
  <img src="png/ex2_bandpass_spectra_core1.png" alt="ex2_bandpass_spectra_core1" style="width:600px;"/>
  <img src="png/ex2_bandpass_spectra_iq.png"    alt="ex2_bandpass_spectra_iq" style="width:600px;"/>
</p>

---
### Example 3: Oversampled second Nyquist zone

In this example, the quad-core ADC oversamples a signal in the second Nyquist zone. Note how the direction of the spectral ramp has changed. This is due to spectral inversion of the even Nyquist zones due to aliasing. The core 1 and IQ spectra look similar (but that does not mean the IQ spectra is correct!).

<p align="center">
  <img src="png/ex3_bandpass_spectra.png"       alt="ex3_bandpass_spectra" style="width:600px;"/>
  <img src="png/ex3_bandpass_spectra_core1.png" alt="ex3_bandpass_spectra_core1" style="width:600px;"/>
  <img src="png/ex3_bandpass_spectra_iq.png"    alt="ex3_bandpass_spectra_iq" style="width:600px;"/>
</p>

---
### Example 4: Oversampled third Nyquist zone

In this example, the quad-core ADC oversamples a signal in the third Nyquist zone. The core 1 and IQ spectra look similar (but that does not mean the IQ spectra is correct!).

<p align="center">
  <img src="png/ex4_bandpass_spectra.png"       alt="ex4_bandpass_spectra" style="width:600px;"/>
  <img src="png/ex4_bandpass_spectra_core1.png" alt="ex4_bandpass_spectra_core1" style="width:600px;"/>
  <img src="png/ex4_bandpass_spectra_iq.png"    alt="ex4_bandpass_spectra_iq" style="width:600px;"/>
</p>

---
### Example 5: Oversampled fourth Nyquist zone

In this example, the quad-core ADC oversamples a signal in the fourth Nyquist zone. The core 1 and IQ spectra look similar (but that does not mean the IQ spectra is correct!).

<p align="center">
  <img src="png/ex5_bandpass_spectra.png"       alt="ex5_bandpass_spectra" style="width:600px;"/>
  <img src="png/ex5_bandpass_spectra_core1.png" alt="ex5_bandpass_spectra_core1" style="width:600px;"/>
  <img src="png/ex5_bandpass_spectra_iq.png"    alt="ex5_bandpass_spectra_iq" style="width:600px;"/>
</p>

---
### Example 6: Oversampled first and second Nyquist zones

In this example, the quad-core ADC oversamples a signal in the first and second Nyquist zones. The "bathtub" spectra appears again because Nyquist zones 1 and 2 alias on top of each other, and Nyquist zone 2 gets frequency reversed. The core 1 and IQ spectra look similar (but that does not mean the IQ spectra is correct!).

<p align="center">
  <img src="png/ex6_bandpass_spectra.png"       alt="ex6_bandpass_spectra" style="width:600px;"/>
  <img src="png/ex6_bandpass_spectra_core1.png" alt="ex6_bandpass_spectra_core1" style="width:600px;"/>
  <img src="png/ex6_bandpass_spectra_iq.png"    alt="ex6_bandpass_spectra_iq" style="width:600px;"/>
</p>

---
### Example 7: Oversampled first and third Nyquist zone

In this example, the quad-core ADC oversamples a signal in the first and third Nyquist zones. The "bathtub" spectra does not appear this time, because Nyquist zones 1 and 3 alias on top of each other, and the spectra in those two Nyquist zones have the same bandpass slope (and frequency sense). The core 1 and IQ spectra look similar (but that does not mean the IQ spectra is correct!).

<p align="center">
  <img src="png/ex7_bandpass_spectra.png"       alt="ex7_bandpass_spectra" style="width:600px;"/>
  <img src="png/ex7_bandpass_spectra_core1.png" alt="ex7_bandpass_spectra_core1" style="width:600px;"/>
  <img src="png/ex7_bandpass_spectra_iq.png"    alt="ex7_bandpass_spectra_iq" style="width:600px;"/>
</p>


