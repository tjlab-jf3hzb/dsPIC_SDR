%
%        b0 + b1*Z^-1 + b2*Z^-2
%  H(z)=-------------------------
%        1  - a1*Z^-1 - a2*Z^-2
%
%
% Coefficients:
%  16bit 1.15format  
%  k: Integer 
%
%
%    -----+-(b0/k)->"+"--------->"+"-->(k)---+---> 
%         |          ^            ^          |
%        Z^-1        |            |         Z^-1
%         |          |            |          |
%         +-(b1/k)->"+"          "+"<-(a1/k)-+ 
%         |          ^            ^          |
%        Z^-1        |            |         Z^-1
%         |          |            |          |
%         +-(b2/k)-->+            +<-(a2/k)--+
%    
%    
 
 
clear;
pkg load signal;

fs=80e3; %[Hz]
fc=3e3; %{Hz]
Rp=2;  %[dB]
Rs=50; %[dB]
ord=4;

bitk=2;   % k=2^bitk


[b,a]=ellip(ord, Rp, Rs, 2*fc/fs);
%[b,a]=cheby1(ord, Rp, 2*fc/fs);


figure(1);
freqz(b,a);

[sos, g]=tf2sos(b,a);
g1=sum(sos(1,1:3))/sum(sos(1,4:6));
g2=sum(sos(2,1:3))/sum(sos(2,4:6));
sos(1,1:3)=sos(1,1:3)./g1;
sos(2,1:3)=sos(2,1:3)./g2;

sos=sos.*( 2^(15-bitk) );
sos=round(sos);

figure(2);
freqz(sos(1,1:3), sos(1,4:6) );
figure(3);
freqz(sos(2,1:3), sos(2,4:6) );


FID=fopen("IIR_LPF.h", "w");
fprintf(FID, "//%dth IIR Elliptic LPF\n" ,ord);
fprintf(FID, "//fs=%d[Hz]\n" ,fs);
fprintf(FID, "//fc=%d[Hz]\n" ,fc);
fprintf(FID, "//Ripple=%d[dB]\n", Rp );
fprintf(FID, "//Att=%d[dB]\n", Rs );


fprintf(FID, "fractional _XDATA(16) IIR_coef0[] =\n" );
fprintf(FID,"{%d,%d,%d,%d,%d,%d};\n", [sos(1,1:3) -1.*sos(1,5:6) bitk] );

fprintf(FID, "fractional _XDATA(16) IIR_coef1[] =\n" );
fprintf(FID,"{%d,%d,%d,%d,%d,%d};\n", [sos(2,1:3) -1.*sos(2,5:6) bitk] );

fclose(FID);

%
% Header file  "IIR_LPF.h"
%
%fractional _XDATA(16) IIR_coef0[] =
%{b0/k, b1/k, b2/k, -a1/k, -a2/k, bitk};
%
%fractional _XDATA(16) IIR_coef1[] =
%{b0/k, b1/k, b2/k, -a1/k, -a2/k, bitk};
%