function [pin] = create_pin(stimdb,onbin,stimpts,irpts,nrep_stim,F0,Fs,t)


pin_t(onbin+1:onbin+stimpts) = sqrt(2)*20e-6*10^((stimdb)/20)*sin(2*pi*F0*t(1:stimpts));
pin_t(onbin+1:onbin+irpts) = pin_t(onbin+1:onbin+irpts).*(0:(irpts-1))/irpts;
pin_t(onbin+(stimpts-irpts):onbin+stimpts)=pin_t(onbin+(stimpts-irpts):onbin+stimpts).*(irpts:-1:0)/irpts;
pin_t(onbin+stimpts:100*Fs*1e-3) = 0;


for pin_k = 1: nrep_stim
    pin((pin_k-1)*100*Fs*1e-3 +1 :(pin_k)*100*Fs*1e-3) = pin_t;
end


end