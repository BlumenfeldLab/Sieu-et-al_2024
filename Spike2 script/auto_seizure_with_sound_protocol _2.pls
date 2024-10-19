            SET    1,1,0           ; Get rate & scaling OK

            VAR    V45,LoopC=0     ; Define variable for section loops
            VAR    V46,RampC=0     ; Define variable for ramp loops
            VAR    V47,DelayC=0    ; Define variable for delay loops
            VAR    V48,Delay2=0    ;  and another one
            VAR    V49,Delay3=0    ;  and another one
            VAR    V50,Delay4=0    ;  and another one
            VAR    V51,Delay5=0    ;  and another one

E0:         DIGOUT [......00]
            DAC    0,0
            DAC    1,0
            DELAY  s(0.996)-1
            HALT                   ; End of this sequence section

EA:     'a  DIGOUT [......00]
            DAC    0,0
            DAC    1,0
            WAVEGO z,T             ; Arm waveform
RA0:        WAVEBR RA0,W           ; Ensure waveform ready
            WAVEST T               ; Trigger waveform
            DELAY  s(0.011)-1
            WAVEST S               ; Stop waveform
            DAC    1,0
            DELAY  s(0.003)-1
            MARK   5               ; Generate digital marker
            DIGOUT [......10]
            DELAY  s(0.047)-1
            DIGOUT [......00]
            DELAY  s(0.003)-1
            DELAY  9999            >Wait 10 to 15 s
            MOVRND DelayC,30
RA1:        ADDI   DelayC,-429496  >Wait 10 to 15 s
            BGT    DelayC,0,RA1    >Wait 10 to 15 s
            DELAY  s(0.02)-1
            JUMP   EA              ; Jump to next section

EB:     'b  DIGOUT [......00]
            DAC    0,0
            DAC    1,0
            DELAY  s(0.006)-1
            DIGOUT [......01]
            DELAY  s(0.009)-1
            DIGOUT [......00]
            DELAY  s(3.479)-1
            JUMP   ED              ; Jump to next section

EC:     'c  DIGOUT [......00]
            DAC    0,0
            DAC    1,0
            DELAY  s(0.036)-1
            DIGOUT [......10]
            DELAY  s(0.959)-1
            HALT                   ; End of this sequence section

ED:     'd  DIGOUT [......00]
            DAC    0,0
            DAC    1,0
            WAVEGO z,T             ; Arm waveform
RD0:        WAVEBR RD0,W           ; Ensure waveform ready
            WAVEST T               ; Trigger waveform
            DELAY  s(0.011)-1
            WAVEST S               ; Stop waveform
            DAC    1,0
            DELAY  s(0.003)-1
            MARK   5               ; Generate digital marker
            DIGOUT [......10]
            DELAY  s(0.047)-1
            DIGOUT [......00]
            DELAY  s(0.003)-1
            DELAY  1999            >Wait 5 to 7 s
            MOVRND DelayC,30
RD1:        ADDI   DelayC,-2147483 >Wait 5 to 7 s
            BGT    DelayC,0,RD1    >Wait 5 to 7 s
            DELAY  s(0.02)-1
            JUMP   ED              ; Jump to next section

ES:     's  DIGOUT [......00]
            DAC    0,0
            DAC    1,0
            DELAY  s(0.007)-1
            DAC    1,0
            DELAY  s(0.099)-1
            DAC    1,0
            DELAY  s(0.888)-1
            HALT                   ; End of this sequence section

