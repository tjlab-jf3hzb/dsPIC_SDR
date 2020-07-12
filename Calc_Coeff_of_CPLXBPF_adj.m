clear;

%-------Modify if you need -------------------------------------------

% Define Pass Band: f1(Hz) to f2(Hz)
fL=200;
fH=2700;

%Sampling frequency
fs=10e3;

%Number of Taps
Ndat=128;

%Set a Window function
%wf=ones(Ndat,1);
%wf=hanning(Ndat);
%wf=hamming(Ndat);
wf=blackman(Ndat);
%wf=bartlett(Ndat);

amp_balance=1.00;
ph_error=0;  %[deg]
delay=-2e-6; %[sec]
%--------------------------------------------------------------------




%---------- Calc. complex FIR coeff. ---------------------
fx=[ [0:Ndat/2-1] [-Ndat/2:-1] ].*(fs/Ndat); 

ka=round( fL/(fs/Ndat) ); kb=round( fH/(fs/Ndat) );
if (ka==0) ka=1; end
if (kb>=(Ndat/2-2)) kb=Ndat/2-2; end

fresH=zeros(1,Ndat/2);
fresH(ka+1:kb+1)=1;

I_fres=[fresH 0 fresH(end:-1:2) ];

ph_v= ( 2.*( (fx>=0)-0.5 ) ).*( -pi/2 + pi*ph_error/180);
ph_v=exp( sqrt(-1).*ph_v );
dl_v=exp( sqrt(-1).*2.*pi.*fx.*delay);

Q_fres=I_fres.*ph_v.*dl_v;

figure(1)
subplot(1,2,1)
plot(fx,abs(I_fres));
subplot(1,2,2)
plot(fx,arg(I_fres));

figure(2)
subplot(1,2,1)
plot(fx,abs(Q_fres));
subplot(1,2,2)
plot(fx,arg(Q_fres));


I_ipr=ifft(I_fres);
Q_ipr=ifft(Q_fres);

Ich=[real(I_ipr(Ndat/2+1:end)) real(I_ipr(1:Ndat/2))];
Qch=[real(Q_ipr(Ndat/2+1:end)) real(Q_ipr(1:Ndat/2))];
Ich=Ich.*wf.'; Qch=Qch.*wf.';


figure(3); plot(Ich); ylabel('FIR coeff. (Real Part)');
figure(4); plot(Qch); ylabel('FIR coeff. (Imaginary Part)');



k=8;
[hr,fr]=freqz(Ich,[1], Ndat*k*2,'whole',fs);
hr=[hr(Ndat*k+1:end); hr(1:Ndat*k)];
fr=fr-fs/2;

figure(5)
subplot(1,2,1)
plot(fr, 20.*log10( abs(hr)) );
xlabel('Frequency(Hz)');
ylabel('Response(dB)');
grid on;
axis([fL-300, fL+300, -80 10]);

subplot(1,2,2)
plot(fr, 20.*log10( abs(hr)) );
xlabel('Frequency(Hz)');
ylabel('Response(dB)');
grid on;
axis([fH-300, fH+300, -80 10]);





%--- Convert format to Q1.15 & Create Header file -------------
Icoeff=round(Ich*32767*0.5);
Qcoeff=round(Qch*32767*0.5*amp_balance);

Icoeff=Icoeff+(Icoeff<0)*65536;
Qcoeff=Qcoeff+(Qcoeff<0)*65536;

FID=fopen("H_filter.h", "w");


fprintf(FID, "#define fs %d\n" , fs);
fprintf(FID, "const fractional H_Re[%d]={\n" , Ndat);

for p=[1:Ndat-1]
fprintf(FID,"0x%04x,\n", Icoeff(p));  
end;
fprintf(FID,"0x%04x\n", Icoeff(Ndat) );
fprintf(FID,"};\n");

fprintf(FID, "const fractional H_Im[%d]={\n" , Ndat);
for p=[1:Ndat-1]
fprintf(FID,"0x%04x,\n", Qcoeff(p) );  
end;
fprintf(FID,"0x%04x\n",  Qcoeff(Ndat) );
fprintf(FID,"};\n");

fclose(FID);

