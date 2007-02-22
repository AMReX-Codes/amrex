; a sample start-up script for a Linux box.  The decomposed=0 is
; essential.  It is important that X be configured in a 24-bit color
; mode -- 16-bit color won't work.

device, retain=2, decomposed=0, true_color=24
print, 'Number of colors allocated is ', !d.n_colors
