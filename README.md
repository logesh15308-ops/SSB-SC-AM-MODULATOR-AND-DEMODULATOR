# SSB-SC-AM-MODULATOR-AND-DEMODULATOR
## AIM

To write a program to perform SSBSC modulation and demodulation using SCI LAB and study its spectral characteristics
## EQUIPMENTS REQUIRED

•	Computer with i3 Processor

•	SCI LAB

Note: Keep all the switch faults in off position

## ALGORITHM
  1.	Define Parameters:
     
    •	Fs: Sampling frequency.
    •	T: Duration of the signal.
    •	Fc: Carrier frequency.
    •	Fm: Frequency of the message signal.
    •	Amplitude: Maximum amplitude of the message signal.
    
  2.	Generate Signals:
     
    •	Message Signal: The baseband signal that will be modulated.
    •	Carrier Signal: A high-frequency signal used for modulation.
    •	Analytic Signal: Constructed using the Hilbert transform to get the in-phase and quadrature components.

  3.	SSBSC Modulation:
     
    •	Modulated Signal: Create the SSBSC signal using the in-phase and quadrature components, modulated by the carrier.

  4.	SSBSC Demodulation:
     
    •	Mixing: Multiply the SSBSC signal with the carrier to retrieve the message signal.
    •	Low-pass Filtering: Apply a low-pass filter to remove high-frequency components and recover the original message signal.

  5.	Visualization:
     
    Plot the message signal, carrier signal, SSBSC modulated signal, and the recovered signal after demodulation.

## PROCEDURE

  •	Refer Algorithms and write code for the experiment.
  
  •	Open SCILAB in System
  
  •	Type your code in New Editor
  
  •	Save the file
   
  •	Execute the code
  
  •	If any Error, correct it in code and execute again
  
  •	Verify the generated waveform using Tabulation and Model Waveform
## MODEL GRAPH
<img width="747" height="338" alt="image" src="https://github.com/user-attachments/assets/dd117c5d-ee32-47c7-946c-df6180b0d33f" />

## PROGRAMfs = 2130 * 10;
Am = 2.6 + 0.1*59;
fm = 213 + 0.1*59;
Ac = 2 * Am;
fc = 45 * 10;

t = 0:1/fs:2/fm;
m = Am*cos(2*%pi*fm*t);
M_f = fft(m);
N = length(M_f);
H = zeros(1, N);

if modulo(N,2)==0 then
    H(1) = 1;
    H(N/2+1) = 1;
    H(2:N/2) = 2;
else
    H(1) = 1;
    H((N+1)/2) = 1;
    H(2:(N+1)/2-1) = 2;
end

mh = real(ifft(M_f .* (-%i*H)));
c = Ac*cos(2*%pi*fc*t);
s = Ac*sin(2*%pi*fc*t);

ssb_usb = m .* c - mh .* s;
ssb_lsb = m .* c + mh .* s;

demod_usb = ssb_usb .* c;
demod_lsb = ssb_lsb .* c;

function y = lowpass(x, cutoff)
    N = length(x);
    X = fft(x);
    f = (0:N-1)*(fs/N);
    H = (f < cutoff);
    y = real(ifft(X .* H));
endfunction

m_rec_usb = lowpass(demod_usb, 2*fm);
m_rec_lsb = lowpass(demod_lsb, 2*fm);

clf;
subplot(3,2,1);
plot(t, ssb_usb);
title("SSB Modulated Signal (USB)");
xlabel("time"); ylabel("amplitude");

subplot(3,2,2);
plot(t, ssb_lsb);
title("SSB Modulated Signal (LSB)");
xlabel("time"); ylabel("amplitude");

subplot(3,2,3);
plot(fftshift(abs(fft(ssb_usb))));
title("Spectrum (USB)");
xlabel("frequency"); ylabel("amplitude");

subplot(3,2,4);
plot(fftshift(abs(fft(ssb_lsb))));
title("Spectrum (LSB)");
xlabel("frequency"); ylabel("amplitude");

subplot(3,2,5);
plot(t, m_rec_usb);
title("Demodulated Message (USB)");
xlabel("time"); ylabel("amplitude");

subplot(3,2,6);
plot(t, m_rec_lsb);
title("Demodulated Message (LSB)");
xlabel("time"); ylabel("amplitude");

## TABULATION![WhatsApp Image 2025-11-20 at 07 07 51_31e33e91](https://github.com/user-attachments/assets/e74d4189-a9d7-46e6-8e82-583008a41a6c)


## OUTPUT![WhatsApp Image 2025-11-20 at 07 08 28_2ca1e105](https://github.com/user-attachments/assets/cb1e0590-0695-40c2-b767-9c5e09b17dc2)


## RESULTThus, the SSB-SC-AM Modulation and Demodulation is experimentally done and the output is verified.


