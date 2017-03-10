maxrmsstr = '5e-4'
maxrms = float(maxrmsstr)
rloudfacstr = '2.7e-2'
rloudfac = float(rloudfacstr)
radsigmafac = 3.
o3sigmafac = 3.
firstsigmafac = 5.
snfac = 3.
snfacstr = str(int(snfac))
lo3leddbr = -4.5
radsnfac = 3.

band = 'C'
if band == 'C':
    freq = 4.9e9
elif band == 'X':
    freq = 8.4e9

niter = 100
linesnfac = 3.
bdecintr = 3.1

setboxdiam = 17 #pixels
numproc = 7
drver = 'dr9'
