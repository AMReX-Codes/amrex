#!/usr/bin/python

import os
import sys
import getopt

def main(argv):

    execname = ""
    inputname = ""
    directory_prefix = ""
    richardsontestloc = ""
    coarsteps = 0
    coarcells = 0
    coarmaxgrid = 0
    trunconly =  ""
    try:
        opts, args = getopt.getopt(argv,"e:i:s:n:m:t:q:H")
    except getopt.GetoptError:
        print 'usage: convtest_serial_driver.py  -e execfile -i inputfile -s coarse_steps -n coarse_cells -m max_grid_cells_coar -t trunconly -q richardsontestloc'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-H':
            print 'usage: convtest_serial_driver.py  -e execfile -i inputfile -s coarse_steps -n coarse_cells -m max_grid_cells_coar -t trunconly  -q richardsontestloc'
            sys.exit()
        elif opt in ("-e"):
            execname = arg
        elif opt in ("-i"):
            inputname = arg
        elif opt in ("-s"):
            coarsteps = int(arg)
        elif opt in ("-n"):
            coarcells  = int(arg)
        elif opt in ("-m"):
            coarmaxgrid  = int(arg)
        elif opt in ("-q"):
            richardsontestloc = arg
        elif opt in ("-t"):
            trunconly  = arg

    print "exec_name is        "   + execname
    print "input_name is   "       + inputname
    print "nx coarsest is      "   + str(coarcells)
    print "max_step coarse is  "   + str(coarsteps)
    print "max_grid coarse is  "   + str(coarmaxgrid)
    print "truncation only is  "   + str(trunconly)
    print "richardson test is  "   + richardsontestloc


    medisteps = 2*coarsteps
    finesteps = 2*medisteps
    medicells = 2*coarcells
    finecells = 2*medicells
    medimaxgrid = 2*coarmaxgrid
    finemaxgrid = 2*medimaxgrid
    if(trunconly == "true"):
        finesteps = 1
        medisteps = 1
        coarsteps = 1

    truncstr = "truncation_error_only=" + trunconly
    ncellcoarstr   = "amr.n_cells=\"" + str(coarcells) + " " + str(coarcells) + " " + str(coarcells) + "\""
    ncellmedistr   = "amr.n_cells=\"" + str(medicells) + " " + str(medicells) + " " + str(medicells) + "\""
    ncellfinestr   = "amr.n_cells=\"" + str(finecells) + " " + str(finecells) + " " + str(finecells) + "\""
    maxgridcoarstr = "amr.maxgrid=\""  + str(coarmaxgrid) + "\""
    maxgridmedistr = "amr.maxgrid=\""  + str(medimaxgrid) + "\""
    maxgridfinestr = "amr.maxgrid=\""  + str(finemaxgrid) + "\""
    maxstepcoarstr = "max_step=\""  + str(coarsteps) + "\""
    maxstepmedistr = "max_step=\""  + str(medisteps) + "\""
    maxstepfinestr = "max_step=\""  + str(finesteps) + "\""


    pltfilecoar = ""
    pltfilemedi = ""
    pltfilefine = ""
    if(coarsteps > 999):
        pltfilecoar = "plt0" + str(coarsteps)
    elif(coarsteps > 99):
        pltfilecoar = "plt00" + str(coarsteps)
    elif(coarsteps > 9):
        pltfilecoar = "plt000" + str(coarsteps)
    else:
        pltfilecoar = "plt0000" + str(coarsteps)

    if(medisteps > 999):
        pltfilemedi = "plt0" + str(medisteps)
    elif(medisteps > 99):
        pltfilemedi = "plt00" + str(medisteps)
    elif(medisteps > 9):
        pltfilemedi = "plt000" + str(medisteps)
    else:
        pltfilemedi = "plt0000" + str(medisteps)

    if(finesteps > 999):
        pltfilefine = "plt0" + str(finesteps)
    elif(finesteps > 99):
        pltfilefine = "plt00" + str(finesteps)
    elif(finesteps > 9):
        pltfilefine = "plt000" + str(finesteps)
    else:
        pltfilefine = "plt0000" + str(finesteps)
    

    print "coar plot file = " + pltfilecoar + ", medi plot file = " + pltfilemedi+ ", fine plot file = " + pltfilefine

    commandcoar = execname + " "  + inputname + " amr.plot_int=100000  " + ncellcoarstr + " " + maxgridcoarstr + " " + maxstepcoarstr + " " + truncstr + " | tee screen.out.coar"
    commandmedi = execname + " "  + inputname + " amr.plot_int=100000  " + ncellmedistr + " " + maxgridmedistr + " " + maxstepmedistr + " " + truncstr + " | tee screen.out.medi"
    commandfine = execname + " "  + inputname + " amr.plot_int=100000  " + ncellfinestr + " " + maxgridfinestr + " " + maxstepfinestr + " " + truncstr + " | tee screen.out.fine"

    cleanupcoar = "mv " + pltfilecoar + " _plt.coar"
    cleanupmedi = "mv " + pltfilemedi + " _plt.medi"
    cleanupfine = "mv " + pltfilefine + " _plt.fine"
    print commandcoar
    os.system(commandcoar)
    print cleanupcoar
    os.system(cleanupcoar)
    print commandmedi
    os.system(commandmedi)
    print cleanupmedi
    os.system(cleanupmedi)
    print commandfine
    os.system(commandfine)
    print cleanupfine
    os.system(cleanupfine)
    
    actualtest = richardsontestloc + " coarFile=_plt.coar mediFile=_plt.medi fineFile=_plt.fine coarError=_err.coar mediError=_err.medi | tee richardson.out"
    print actualtest
    os.system(actualtest)

if __name__ == "__main__":
   main(sys.argv[1:])

