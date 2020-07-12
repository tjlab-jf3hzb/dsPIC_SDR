/* 
 * File:   main.c
 * 
 * Author: tj. Uebo/ JF3HZB
 * 
 * " SDR  for SSB"
 * 
 *    July, 09, 2020
 *
 *   XC16 v1.40  with "libdsp-elf.a"
 *   MPLAB 5.25
 */

#include <xc.h>
#include <dsp.h>

// DSPIC33FJ64GP802 Configuration Bit Settings
// FBS
#pragma config BWRP = WRPROTECT_OFF     // Boot Segment Write Protect (Boot Segment may be written)
#pragma config BSS = NO_FLASH           // Boot Segment Program Flash Code Protection (No Boot program Flash segment)
#pragma config RBS = NO_RAM             // Boot Segment RAM Protection (No Boot RAM)
// FSS
#pragma config SWRP = WRPROTECT_OFF     // Secure Segment Program Write Protect (Secure segment may be written)
#pragma config SSS = NO_FLASH           // Secure Segment Program Flash Code Protection (No Secure Segment)
#pragma config RSS = NO_RAM             // Secure Segment Data RAM Protection (No Secure RAM)
// FGS
#pragma config GWRP = OFF               // General Code Segment Write Protect (User program memory is not write-protected)
#pragma config GSS = HIGH               // General Segment Code Protection (High Security Code Protection is Enabled)
// FOSCSEL
#pragma config FNOSC = FRCPLL           // Oscillator Mode (Internal Fast RC (FRC) w/ PLL)
#pragma config IESO = ON                // Internal External Switch Over Mode
// FOSC
#pragma config POSCMD = NONE            // Primary Oscillator Source (Primary Oscillator Disabled)
#pragma config OSCIOFNC = ON            // OSC2 Pin Function (OSC2 pin has digital I/O function)
#pragma config IOL1WAY = ON             // Peripheral Pin Select Configuration (Allow Only One Re-configuration)
#pragma config FCKSM = CSDCMD           // Clock Switching and Monitor (Both Clock Switching and Fail-Safe Clock Monitor are disabled)
// FWDT
#pragma config WDTPOST = PS32768        // Watchdog Timer Postscaler (1:32,768)
#pragma config WDTPRE = PR128           // WDT Prescaler (1:128)
#pragma config WINDIS = OFF             // Watchdog Timer Window (Watchdog Timer in Non-Window mode)
#pragma config FWDTEN = OFF             // Watchdog Timer Enable (Watchdog timer enabled/disabled by user software)
// FPOR
#pragma config FPWRT = PWR128           // POR Timer Value (128ms)
#pragma config ALTI2C = OFF             // Alternate I2C  pins (I2C mapped to SDA1/SCL1 pins)
// FICD
#pragma config ICS = PGD2               // Comm Channel Select (Communicate on PGC2/EMUC2 and PGD2/EMUD2)
#pragma config JTAGEN = OFF             // JTAG Port Enable (JTAG is Disabled)

//#define LED  LATAbits.LATA2
#define UL    PORTBbits.RB7
#define LSB 1
#define USB 0
volatile unsigned char f_SIDEBAND=LSB;



fractional _YBSS(16) Ich=0, Qch=0;
fractional _YBSS(16) Isig=0, Qsig=0, Iout=0, Qout=0;
fractional _YBSS(16) Re_in=0, Im_in=0, Re_out=0,Im_out=0;
fractional V0, V1;
char f_ANch=0;


#include "H_filter.h"
//---Complex  FIR-------------------
#define Ndat 128
#define Bdat 7
fractional _XDATA(2*Ndat) Re_coeff[Ndat];
FIRStruct	CPLXFilter_Re;
fractional _YBSS(2*Ndat) delay_Re[Ndat];

fractional _XDATA(2*Ndat) Im_coeff[Ndat];
FIRStruct	CPLXFilter_Im;
fractional _YBSS(2*Ndat) delay_Im[Ndat];


//--- IIR -----------------------------------
#include "IIR_LPF.h"
fractional _YBSS(16) Z_Re_0[5];
fractional _YBSS(16) Z_Re_1[5];
fractional _YBSS(16) Z_Re_2[5];
fractional _YBSS(16) Z_Re_3[5];
fractional _YBSS(16) Z_Im_0[5];
fractional _YBSS(16) Z_Im_1[5];
fractional _YBSS(16) Z_Im_2[5];
fractional _YBSS(16) Z_Im_3[5];


/*------------------------------------------------------------------------------------------
 * -----------------------------------------------------------------------------------------
 *           main
 * -----------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------*/
int main(int argc, char** argv) {
  
	//------- Clock -------------------------------------
	CLKDIVbits.PLLPRE = 0;		//N1=1/2
	CLKDIVbits.PLLPOST = 0;		//N2=1/2
	//PLLFBD = 42;				//M=42+2, Fin=M*7.37(FSC)/2/2/2 = 40.535MHz(Fcy)
    PLLFBD = 50;				//M=50+2, Fin=M*7.37(FSC)/2/2/2 = 47.9MHz(Fcy)
    //-- Setup CPLX FIR Filter ---------------------------------------------------------
    FIRStructInit( &CPLXFilter_Re, Ndat, (fractional *)Re_coeff, COEFFS_IN_DATA, delay_Re );
	FIRDelayInit( &CPLXFilter_Re );   
    FIRStructInit( &CPLXFilter_Im, Ndat, (fractional *)Im_coeff, COEFFS_IN_DATA, delay_Im );
	FIRDelayInit( &CPLXFilter_Im );

   
	//------ Port ----------------------------------
	//TRISAbits.TRISA2 = 0;
    //LED=0;
	//LATB  = 0x0000;
	TRISB = 0xFFFF;

    //--- TIMER1  --------------------------------------------   
    //#define Ndiv_fsi 480        //for fs=100kHz
    #define Ndiv_fsi 600        //for fs=80kHz
	IEC0bits.T1IE = 1;			//enable TIMER1 INT	
	T1CONbits.TCKPS = 0;		//prescaler : 1/1
	PR1 = (Ndiv_fsi-1);         // Fcy/(Ndiv_fsi)
	TMR1 = 0;	
    IPC0bits.T1IP=6;	
    //T1CONbits.TON =1;			//start TIMER1
       
	//--- ADC ------------------------------------------------
	AD1CON1 = 0x07E0;		// 12bit Signed Fractional, Manual Sample & Auto Start Conv.  
	AD1CON2 = 0;            // Vref: AVdd, AVss    
	AD1CON3 = 0x0802;		// Tsmp=8*Tad, Tad = 3*Tcy,  Tconv=(8+14)*Tad= 66/Fcy   
    AD1CHS0=0x0404;   
	AD1PCFGL = 0xFFCC;		// Use AN0,AN1,AN4,AN5   
	AD1CON1bits.ADON = 1;	// Enable AD
    f_ANch=0;           
    IFS0bits.AD1IF=0;
    IPC3bits.AD1IP=7;
    IEC0bits.AD1IE=1;
    //AD1CON1bits.SAMP=1;   // Start sample
        
    //--- DAC ------------------------------------------------
    ACLKCON=0x0700;  // ACLK=Fcy/1
    
    //---- fs -----------------------------------
    #define Ndiv 75  // Fcy/64/1/Ndiv=10kHz
    //#define Ndiv 60  // Fcy/64/1/Ndiv=12.5kHz

    DAC1STATbits.ROEN = 1; // Enable Right Channel DAC Output
    DAC1STATbits.LOEN = 1; // Enable Left Channel DAC Output 
    DAC1STATbits.RITYPE = 0; // Interrupt Right Channel if FIFO is not Full
    DAC1STATbits.LITYPE = 0; // Interrupt Left Channel if FIFO is not Full 
    DAC1CONbits.DACFDIV = Ndiv-1;  // Divide Clock by Ndiv
    DAC1CONbits.FORM = 1; // Data Format is Signed
    DAC1DFLT = 0x0000; // Default value when the FIFO is empty
    DAC1LDAT = 0x0000;
    DAC1RDAT = 0x0000;   
    DAC1CONbits.DACEN = 1; // Enable DAC1 Module
    IFS4bits.DAC1LIF = 0; // Clear Left Channel Interrupt Flag  
    IFS4bits.DAC1RIF = 0; // Clear Right Channel Interrupt Flag
    IPC19bits.DAC1LIP = 5;
    IPC19bits.DAC1RIP = 0;

    IEC4bits.DAC1LIE = 1; // Disable Left Channel Interrupt
    IEC4bits.DAC1RIE = 0; // Disable Right Channel Interrupt
    //DAC1LDAT=0;
    //DAC1RDAT=0;
    T1CONbits.TON =1;			//start TIMER1

    int k;
    for (k = 0; k < Ndat; k++) {
        Re_coeff[k] = H_Re[k];
        Im_coeff[k] = H_Im[k];
    }

    
    /*--------------------------------------------------------------------------
          Main Loop
     -------------------------------------------------------------------------*/
    while (1) {
        f_SIDEBAND=UL;
    }
    /*--------------------------------------------------------------------------
          End of Main Loop
     -------------------------------------------------------------------------*/  
    return (EXIT_SUCCESS);
}







void __attribute__((interrupt, no_auto_psv))_ADC1Interrupt(void){

    IFS0bits.AD1IF=0;  
    if(f_ANch==0){
        Ich = ADC1BUF0;
        AD1CHS0=0x0505;
        AD1CON1bits.SAMP=1; // Start sampling
    }
    else if(f_ANch==1)
    {
        Qch = ADC1BUF0;
        AD1CHS0=0x0000;
        AD1CON1bits.SAMP=1; // Start sampling
    }
    else if(f_ANch==2){
        V0 = ADC1BUF0;   // not used
        AD1CHS0=0x0101;
        AD1CON1bits.SAMP=1; // Start sampling
    }    
    else if(f_ANch==3){
        V1 = ADC1BUF0;   // not used
        AD1CHS0=0x0404;
    }        
    f_ANch++; f_ANch&=3;    
}




void __attribute__((interrupt, no_auto_psv))_T1Interrupt(void)
{
    fractional y;
    Isig = Ich;
    Qsig = Qch;
    AD1CON1bits.SAMP = 1; // Start sampling
    
/*----- Decimation Filter -----------------------------
%
%        b0 + b1*Z^-1 + b2*Z^-2
%  H(z)=-------------------------
%        1  + a1*Z^-1 + a2*Z^-2
%
%
% Coefficients:
%  16bit 1.15format  
%  k: Integer (k=4)
%
%
%     x --+--(b0/k)-->"+"--->+"--->(k)---+---> y
%         |            ^     ^           |
%        Z^-1          |     |          Z^-1
%         |            |     |           |
%         +--(b1/k)-->"+"   "+"--(a1/k)--+ 
%         |            ^     ^           |
%       Z^-1           |     |          Z^-1
%         |            |     |           |
%         +--(b2/k)-->"+"    + --(a2/k)--+
%    
%
*/
     //--- I input --------------------------------------------
    //-- 1st ------------------------------------------
    Z_Re_0[0]=Isig;
    y=VectorDotProduct(5, &Z_Re_0[0], &IIR_coef0[0]);
    y<<=2;
    Z_Re_0[2]=Z_Re_0[1];
    Z_Re_0[1]=Z_Re_0[0];
    Z_Re_0[4]=Z_Re_0[3];
    Z_Re_0[3]=y;
        
    //-- 2nd ------------------------------------------
    Z_Re_1[0]=y;
    y=VectorDotProduct(5, &Z_Re_1[0], &IIR_coef1[0]);
    y<<=2;
    Z_Re_1[2]=Z_Re_1[1];
    Z_Re_1[1]=Z_Re_1[0];
    Z_Re_1[4]=Z_Re_1[3];
    Z_Re_1[3]=y;
   
    Iout=y;
    
    
    //--- Q input --------------------------------------------    
    //-- 1st ------------------------------------------
    Z_Im_0[0]=Qsig;
    y=VectorDotProduct(5, &Z_Im_0[0], &IIR_coef0[0]);
    y<<=2;
    Z_Im_0[2]=Z_Im_0[1];
    Z_Im_0[1]=Z_Im_0[0];
    Z_Im_0[4]=Z_Im_0[3];
    Z_Im_0[3]=y;
         
    //-- 2nd ------------------------------------------
    Z_Im_1[0]=y;
    y=VectorDotProduct(5, &Z_Im_1[0], &IIR_coef1[0]);
    y<<=2;
    Z_Im_1[2]=Z_Im_1[1];
    Z_Im_1[1]=Z_Im_1[0];
    Z_Im_1[4]=Z_Im_1[3];
    Z_Im_1[3]=y;
    
    Qout=y;

    IFS0bits.T1IF=0;
}






/*------------------------------------------------------------------------------
    Demodulation
------------------------------------------------------------------------------*/
int STagc=0;
long int AGC=0;
long int C_AGC=0;
void __attribute__((interrupt, no_auto_psv))_DAC1LInterrupt(void)
{
    int af;     
    IFS4bits.DAC1LIF = 0;  // Clear Left Channel Interrupt Flag

    Re_in = Iout;
    Im_in = Qout;
    //----- ComplexFIR Filter  ----------------      
    FIR(1, &Re_out, &Re_in, &CPLXFilter_Re);
    FIR(1, &Im_out, &Im_in, &CPLXFilter_Im);

    if (f_SIDEBAND == LSB) {
        af = Re_out - Im_out; //LSB
    } else {
        af = Re_out + Im_out; //USB
    }
            
    //AGC-----------------------------------------------------
    long afl=(long)af;
    long gain=60000 - (AGC>>13);
    
    afl=( afl * gain )>>13;
    
    long aff=afl;
    if( aff>=32767) aff=32767;
    if( aff<=-32767) aff=-32767;      
    DAC1RDAT=(int)(aff);  

    //signal level
    long sig_level;
    if(afl<0) sig_level=(long)(-afl);
        else  sig_level=(long)( afl);
    
    long int datagc; 
    datagc=sig_level<<15;   


    switch (STagc) {
        
        case 0: //State of Attack
            if(AGC<datagc){
                C_AGC+=(datagc-C_AGC)>>4;
                STagc=0;
            }else
            {
                STagc=1;
            }
	        AGC=C_AGC + ((datagc-C_AGC)>>7);
            break;
            
        case 1: //State of Release
            if(AGC>datagc){
                C_AGC-=40000;
                if(C_AGC<0) C_AGC=0;
                STagc=1;
            }else
            {
                STagc=0;
            }
			AGC=C_AGC;
            break;

        default:
            STagc = 0;
    }   
    

    //-----------------------------------------------------------
    
    DAC1LDAT=0;   
    
}

/*
void __attribute__((interrupt, no_auto_psv))_DAC1RInterrupt(void)
{
    IFS4bits.DAC1RIF = 0;
}
*/
