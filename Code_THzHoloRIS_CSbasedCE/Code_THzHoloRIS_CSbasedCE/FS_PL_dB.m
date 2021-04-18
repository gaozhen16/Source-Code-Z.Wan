function PL_dB = FS_PL_dB(f,d)
    c = 3e8;
    PL_dB = -20*log10(c/4/pi/f/d);
end