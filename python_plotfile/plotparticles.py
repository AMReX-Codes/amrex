#!/usr/bin/env python
#
# routine to read particle data and create various plots and/or animations.
# 
# 2012-3-01 R. Orvedahl
#

import math
import numpy
import sys
import parseparticles
import pylab
import os
import random
import string
import time as tm
import plotsinglevar_parts

#-----------------------------------------------------------------------------
def main(files):

    # get start time for timing the program
    starttime = tm.time()

    # to allow for an inputs file of variables
    # dict to hold all variables
    variables = {}

    # to run the script in the background set flag to true and use:
    #    nohup plotparticles.py timestamp_* &
    variables['run_in_background'] = True

    # use plotsinglevar_parts.py to plot particles on top of another
    # variable:
    variables['component_plots'] = True

    # 1st and 2nd components to plot using plotsinglevar_parts.py:
    variables['component_one'] = "tfromp"
    variables['component_two'] = ""

    # what kind of movie to make mp4 or avi:
    variables['movie'] = "mp4"

    # save plots as an eps:
    variables['eps'] = False

    # compute turnover time of particles:
    variables['turnover_time'] = False

    # number of grids to use in calculating the turnover time
    variables['num_grids'] = "3"

    # make animations in time
    variables['animation'] = False

    # graph positions
    variables['position_plots'] = False

    # what to graph
    variables['y_vs_x'] = True
    variables['y_vs_t'] = False
    variables['x_vs_t'] = False

    # what fraction of the particles to plot:
    # program will plot frac*(num particles)
    variables['frac'] = "0.1"

    # if not running in background ask whether or not 
    #  to read variable values from an inputs file
    if variables['run_in_background']:
        read_input_file = 'n'    
    else:
        read_input_file = raw_input("\nRead Variables From Input File?:  ")

    if (read_input_file[0] == 'y' or read_input_file[0] == 'Y'):

        filename = raw_input("\n   Input Filename:  ")

        # read variables from input file
        read_from_input_file(variables,filename)

    # maximum value of frac is 1.0
    frac = min(1.0,float(variables['frac']))

    print "\nCounting Particles:"

    # this returns a dict whose keys are a unique identifier (based on 
    # id and CPU) and values are the actual particle objects
    particlesDict = parseparticles.parseParticleFile(files)

    # get just the particle objects
    particles = particlesDict.values()

    print "   Unique Particles:", len(particles)
    print "   Only Using:", int(frac*len(particles))

    # output directory for all images
    outputDir = "Images"

    # current directory
    topDir = os.getcwd()

    # create outputDir if it does not exist already
    if (not os.path.isdir(outputDir)):
        os.mkdir(outputDir)

    # clear plot
    pylab.clf()

    # global domain extrema
    xmin = 1.e33
    xmax = -1e33
    ymin = 1.e33
    ymax = -1.e33

    # Find xmin,xmax,ymin,ymax for all particles and
    # take into account the periodic boundary conditions in x direction
    # this causes lines when the particle jumps from one edge of the 
    # domain to the other. loop through the particles and find where 
    # the x jumps a significant amount and split the x plots according to 
    # these positions.
    x_jump_indices = {}
    n=0
    while (n < int(frac*len(particles))):

        # get numpy arrays containing the time and coordinate
        # information for particle n
        coords, time = particles[n].flatten()

        # find the extrema of the particles over time to make all frames the
        # same size
        xmin = min(xmin,min(coords[0,:]))
        xmax = max(xmax,max(coords[0,:]))
        ymin = min(ymin,min(coords[1,:]))
        ymax = max(ymax,max(coords[1,:]))

        # reset indices list
        indices = []

        i=0
        while (i < len(coords[0,:])-1):

            if (abs(coords[0,i] - coords[0,i+1]) >= 1.0e7):

                # holds the index where a jump across the boundary occurs
                indices.append(i)

            i += 1

        # each particle has a list of indices where 
        # it crosses the boundary
        x_jump_indices[n] = indices

        n += 1

    # COMPONENT PLOTS: plot particle data on top of another variable or
    # on top of two variables using the plotsinglevar_plots.py routine
    if (variables['component_plots']):

       print "Plotting Particles and Components:"
       print "\n   Plotfiles from ./plt*"

       # IDEA:
       # -loop over all pltfiles and get the time of each plotfile
       # -plot the component data from the pltfile
       # -overplot the particle data from that time

       # list to hold all the pltfile directories
       pltFiles = []
       
       # find all directories with name plt?????
       for (path,dirs,files) in os.walk(topDir):
          for directory in dirs:
             index = string.find(directory, "plt")
             if not(index < 0):
                pltFiles.append(directory)

       # dictionary to hold time of each plot file
       pltFile_times = {}
       
       # get time for each plot file
       for plotfile in pltFiles:

          time_pltfile = []

          # use fsnapshot routine to find time of plotfile
          plot_file_time = \
               plotsinglevar_parts.fsnapshot.fplotfile_get_time(plotfile)

          # use a list so it can hold both the time and the index of that time
          time_pltfile.append(plot_file_time)

          pltFile_times[plotfile] = time_pltfile

       # return time and coord data for particle 0
       coords, time = particles[0].flatten()

       # the time from particles[0].flatten() do not exactly match (to 
       # machine precision) the times from the fsnapshot fortran routine.
       # use a tolerance to find the index associated to the time
       time_tol = 1.0e-4

       nstep = 0
       while (nstep < len(time[:])):

           # loop over all plotfiles
           for key, value in sorted(pltFile_times.iteritems()):

              # compare the times
              if (abs(time[nstep] - pltFile_times[key][0]) < time_tol):

                 # attach the time index to the dictionary
                 pltFile_times[key].append(nstep)

           nstep += 1

       # component data and particle data will be plotted using 
       # plotsinglevar_parts.py (which uses fsnapshot)
       count = 1
       for plotfile, value in sorted(pltFile_times.iteritems()):

           pylab.clf()

           # generate a specific output name for the image
           if (variables['component_two'] != ""):
               outFile = variables['component_one'] + str("_") + \
                         variables['component_two'] + str("_particles_") + \
                         plotfile
           else:
               outFile = variables['component_one'] + str("_particles_") + \
                         plotfile

           if (variables['eps']):
               outFile += ".eps"
               eps_send = 1
           else:
               outFile += ".png"
               eps_send = 0

           # plot component data
           # calling sequence: plotfile, comp1, comp2, outfile, log,
           #                   minval, maxval, minval2, maxval2, eps, 
           #                   dpi, origin, annotation, particles, time_ind
           plotsinglevar_parts.do_plot(plotfile, 
                 variables['component_one'], variables['component_two'],
                 outputDir+"/"+outFile, 0, None, None, None, None, 
                 eps_send, 100, 0, "", particles, 
                 pltFile_times[plotfile][1])

           if (count % 5 == 0):
              print "      plotfile: "+str(count)+\
                    " of "+str(len(pltFile_times))

           count += 1

       #
       # Animate the particle positions:
       #
       # output file to contain list of all image files
       list_name = "list.txt"
       if os.path.isfile(list_name):
          os.remove(list_name)

       if (variables['eps']):
           os.system("ls -v "+outputDir+"/"+\
              variables['component_one']+"*_particles_*.eps > "+list_name)
       else:
           os.system("ls -v "+outputDir+"/"+\
              variables['component_one']+"*_particles_*.png > "+list_name)

       # use mencoder to generate an animation from the list
       fps = "15" # frames per sec

       # mencoder cannot handle eps files
       if not(variables['eps']):

          if (variables['movie']=="mp4"):

             movie_output3 = "Component-Animation.mp4"

             if os.path.isfile(movie_output3):
                os.remove(movie_output3)

             string1 = "mencoder mf://@"+list_name+" -o " + movie_output3
             string2 = " -of lavf -lavfopts format=mp4 -ss 1 -ovc "
             string3 = "x264 -x264encopts "
             string4 = "crf=20.0:nocabac:level_idc=30:"
             string5 = "global_header:threads=2 -fps " + fps

             os.system(string1+string2+string3+string4+string5)

          elif (variables['movie']=="avi"):

             movie_output3 = "Component-Animation.avi"

             if os.path.isfile(movie_output3):
                os.remove(movie_output3)

             string1 = "mencoder mf://@"+list_name+" -o " + movie_output3
             string2 = " -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=3000"
             string3 = ":vhq -mf type=png -fps " +fps

             os.system(string1+string2+string3)

          else:
             print "\n---ERROR: unknown movie format---\n"

          os.remove(list_name)

       else:
          print "\n---ERROR: mencoder cannot animate eps images---"
          print   "          list of image files: " + list_name +"\n"


    # TURNOVER TIME for all particles
    if (variables['turnover_time']):

       print "Calculate Turnover Time:"
       print "   Only Using: " + str(len(particles))

       # divide domain into how many grids
       num_grids = int(variables['num_grids'])

       # complete cycle up means go from y_bot to y_top
       y_top = float(num_grids-1)/float(num_grids)*(ymax-ymin) + ymin
       y_bot = (ymax-ymin)/float(num_grids) + ymin

       # dictionaries/lists to hold data
       turnover_up_times = {}    # up data of each particle
       turnover_down_times = {}  # down data of each particle
       counts = {}
       tot_up_count = 0
       tot_down_count = 0
       # put all up data in one list and all down data in one list
       up_data_all_particles = []
       down_data_all_particles = []

       # loop over each particle
       n = 0
       while (n < len(particles)):

           # coordinate and time info for nth particle
           coords, time = particles[n].flatten()

           # reset various counters and lists
           time_in_bot = 0
           time_in_top = 0
           flag_bot = 0
           flag_top = 0
           up_count = 0
           down_count = 0
           count = []
           up_period = []
           down_period = []

           # sometimes len(coords[1,:]) is less 
           # than len(particles[0].history)
           if (len(coords[1,:]) == len(particles[0].history)):
              loop_max = len(particles[0].history)
           else:
              loop_max = len(coords[1,:])

           # loop over time
           nstep = 0
           while (nstep < loop_max):

               # particle is below y_bot
               if (coords[1,nstep] < y_bot):

                   flag_bot = 1
                   time_in_bot = nstep

               # particle is above y_top
               if (coords[1,nstep] > y_top):

                   flag_top = 1
                   time_in_top = nstep

               # particle has made journey either
               # top --> bottom or bottom --> top
               if ((flag_bot == 1) and (flag_top == 1)):

                   # particle went down
                   if (time_in_top < time_in_bot):
                      down = time[time_in_bot]-time[time_in_top]
                      down_period.append(down)
                      down_data_all_particles.append(down)
                      flag_top = 0
                      down_count += 1
                      tot_down_count += 1

                   # particle went up
                   else:
                      up = time[time_in_top]-time[time_in_bot]
                      up_period.append(up)
                      up_data_all_particles.append(up)
                      flag_bot = 0
                      up_count += 1
                      tot_up_count += 1

               nstep += 1

           # debugging:
           count.append(up_count)
           count.append(down_count)
           counts[n] = count

           turnover_up_times[n] = up_period
           turnover_down_times[n] = down_period

           n += 1

       # sanity check/debugging purposes:
#       n=0
#       while n < len(particles):
#           if (len(turnover_up_times[n][:]) != counts[n][0]):
#               print "---ERROR: Up counts differ for particle "+str(n)
#           if (len(turnover_down_times[n][:]) != counts[n][1]):
#               print "---ERROR: Down counts differ for particle "+str(n)
#           n += 1

       #--------------------------------------------------------
       #  Generate Histograms
       #--------------------------------------------------------
       print "Generating Histograms:"

       average_up_time = sum(up_data_all_particles)/ \
                         float(len(up_data_all_particles))

       average_down_time = sum(down_data_all_particles)/ \
                         float(len(down_data_all_particles))

       turnover_time = average_up_time + average_down_time

       # individual histograms
       # (one histogram per particle)
       n=0
       while (n < len(particles)):

           pylab.clf()

           pylab.hist([turnover_up_times[n][:],turnover_down_times[n][:]],
                      bins=50,range=[0,250],color=['b','r'],
                      label=['Up','Down'],histtype='bar')

           pylab.legend()
           pylab.ylabel('Number')
           pylab.xlabel('Bins of 5 sec')
           pylab.figtext(0.65,0.73,'Number of Grids: '+str(num_grids))

           pylab.title('Histogram - Up/Down Time for Particle ' + str(n))

           outFile = outputDir+"/hist_part_%04d_grids_" % n+str(num_grids)

           if (variables['eps']):
               outFile += ".eps"
           else:
               outFile += ".png"
 
           pylab.savefig(outFile)

           n += 1

       # all particles on one histogram
       pylab.clf()
       pylab.hist([up_data_all_particles,down_data_all_particles], 
                   bins=50,range=[0,250],color=['b','r'],
                   label=['Up','Down'],histtype='bar')

       pylab.legend()
       pylab.ylabel('Number')
       pylab.xlabel('Bins of 5 sec')
       pylab.figtext(0.65,0.73,'Number of Grids: '+str(num_grids))
       pylab.figtext(0.6,0.69,"Avg Turnover: "+
                  str(int(1000*turnover_time)/1000.)+" sec")

       pylab.title('Histogram - Up/Down Time for All Particles')

       outFile = outputDir+"/histogram_grids_"+str(num_grids)

       if (variables['eps']):
           outFile += ".eps"
       else:
           outFile += ".png"

       pylab.savefig(outFile)

       print
       print "   Histogram Results:"
       print "      Number of Grids: "+str(num_grids)
       print "      Avg Up Time  : "+str(average_up_time)+" sec"
       print "      Avg Down Time: "+str(average_down_time)+" sec"
       print
       print "      Avg Turnover Time: "+str(turnover_time)+" sec"

    
    #--------------------------------------------------------
    # ANIMATION of positions
    #
    # These animations will plot the entire path of a particle each frame.
    #  Every frame, an entire path is added to the picture. The first
    #  frame has one particles path, the second frame has two particle
    #  paths and the last frame has all particle paths.
    #
    #--------------------------------------------------------
    if (variables['position_plots']):

       # give value to nstep inorder to call plotting routine
       # it has no other use in this loop
       nstep = 2

       print "Creating Images:"

       # plot x vs. y
       if (variables['y_vs_x']):

           print " plotting y vs x"
           n = 0
           while (n < int(frac*len(particles))):

               # get numpy arrays containing the time and coordinate
               # information for particle 0
               coords, time = particles[n].flatten()

               pylab.scatter([coords[0,0]], [coords[1,0]], marker="x")

               # call plotting routine
               plotting(variables, x_jump_indices, n, 
                           nstep, time, coords)

               pylab.xlabel("x")
               pylab.ylabel("y")
               pylab.axis([xmin,xmax,ymin,ymax])

               outFile = "particle_paths_%04d" % n

               if (variables['eps']):
                   outFile += ".eps"
               else:
                   outFile += ".png"
  
               if (os.path.isfile(outputDir+"/"+outFile)):
                   os.remove(outputDir+"/"+outFile)

               pylab.savefig(outputDir+"/"+outFile)

               n += 1

       # plot y vs time
       elif (variables['y_vs_t']):

           print " plotting y vs t"
           n = 0
           while (n < int(frac*len(particles))):

               # get numpy arrays containing the time and coordinate
               # information for particle 0
               coords, time = particles[n].flatten()

               pylab.scatter([time[0]], [coords[1,0]], marker="x")

               # call plotting routine
               plotting(variables, x_jump_indices, n, 
                           nstep, time, coords)

               pylab.xlabel("time")
               pylab.ylabel("y")
               pylab.axis([min(time[:]),max(time[:]),ymin,ymax])

               outFile = "particle_paths_%04d" % n

               if (variables['eps']):
                   outFile += ".eps"
               else:
                   outFile += ".png"

               if (os.path.isfile(outputDir+"/"+outFile)):
                   os.remove(outputDir+"/"+outFile)

               pylab.savefig(outputDir+"/"+outFile)

               n += 1

       # plot x vs time
       elif (variables['x_vs_t']):

           print " plotting x vs t"
           n = 0
           while (n < int(frac*len(particles))):

               # get numpy arrays containing the time and coordinate
               # information for particle 0
               coords, time = particles[n].flatten()

               pylab.scatter([time[0]], [coords[0,0]], marker="x")

               # call plotting routine
               plotting(variables, x_jump_indices, n, 
                           nstep, time, coords)
   
               pylab.xlabel("time")
               pylab.ylabel("x")
               pylab.axis([min(time[:]),max(time[:]),xmin,xmax])

               output = "particle_paths_%04d" % n

               if (variables['eps']):
                   output += ".eps"
               else:
                   output += ".png"

               if (os.path.isfile(outputDir+"/"+output)):
                   os.remove(outputDir+"/"+output)

               pylab.savefig(outputDir+"/"+output)

               n += 1


       # output a file containing all images
       list_name = "list.txt"
       if os.path.isfile(list_name):
          os.remove(list_name)

       if (variables['eps']):
           os.system("ls -v "+outputDir+"/particle_paths_*.eps > "+list_name)
       else:
           os.system("ls -v "+outputDir+"/particle_paths_*.png > "+list_name)

       # use mencoder to generate an animation from the list
       fps = "5" # frames per sec

       # mencoder cannot handle eps files
       if not(variables['eps']):
 
          if (variables['movie']=="mp4"):

             movie_output2 = "Position-Animation.mp4"

             if os.path.isfile(movie_output2):
                os.remove(movie_output2)

             string1 = "mencoder mf://@"+list_name+" -o " + movie_output2
             string2 = " -of lavf -lavfopts format=mp4 -ss 1 -ovc "
             string3 = "x264 -x264encopts "
             string4 = "crf=20.0:nocabac:level_idc=30:"
             string5 = "global_header:threads=2 -fps " + fps

             os.system(string1+string2+string3+string4+string5)

          elif (variables['movie']=="avi"):

             movie_output2 = "Position-Animation.avi"

             if os.path.isfile(movie_output2):
                os.remove(movie_output2)

             string1 = "mencoder mf://@"+list_name+" -o " + movie_output2
             string2 = " -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=3000"
             string3 = ":vhq -mf type=png -fps " +fps

             os.system(string1+string2+string3)

          else:
             print "\n---ERROR: unknown movie format---\n"

          os.remove(list_name)

       else:
          print "\n---ERROR: mencoder cannot animate eps images---"
          print   "          list of image files: " + list_name +"\n"

    # ANIMATION in time
    if (variables['animation']):

       print "Creating Animation:"

       if (variables['y_vs_x']):
          print " animating y vs x"

       elif (variables['y_vs_t']):
          print " animating y vs t"

       elif (variables['x_vs_t']):
          print " animating x vs t"
       
       # make an animation -- note: this assumes that all particles exist
       # at all timesteps
       print "  Total timesteps: "+str(len(particles[0].history))

       nstep = 0
       while (nstep < len(particles[0].history)):
        
           pylab.clf()
        
           n = 0
           while (n < int(frac*len(particles))):

               # get numpy arrays containing the time and coordinate
               # information for particle n
               coords, time = particles[n].flatten()

               # label particle position with an "x"
               pylab.scatter([particles[n].history[nstep].xyz[0]],
                          [particles[n].history[nstep].xyz[1]] ,
                          marker="x")

               # call plotting routine
               plotting(variables, x_jump_indices, n, 
                           nstep, time, coords)

               n += 1

           # axis labels
           if (variables['y_vs_x']):
              pylab.xlabel("x")
              pylab.ylabel("y")
              pylab.axis([xmin,xmax,ymin,ymax])

           elif (variables['x_vs_t']):
              pylab.xlabel("time")
              pylab.ylabel("x")
              pylab.axis([min(time[0:]),max(time[0:]),xmin,xmax])

           elif (variables['y_vs_t']):
              pylab.xlabel("time")
              pylab.ylabel("y")
              pylab.axis([min(time[0:]),max(time[0:]),ymin,ymax])


           output = "particles_%04d" % nstep

           if (variables['eps']):
               output += ".eps"
           else:
               output += ".png"

           if (os.path.isfile(outputDir+"/"+output)):
               os.remove(outputDir+"/"+output)

           # one image per timestep
           pylab.savefig(outputDir+"/"+output)
           
           # update time:
           if nstep % 500 == 0:
              print "   nstep = "+str(nstep)
           nstep += 100

       if (variables['eps']):
           os.system("ls -v "+outputDir+"/particles_*.eps > "+list_name)
       else:
           os.system("ls -v "+outputDir+"/particles_*.png > "+list_name)

       # output a list of particle files
       list_name = "list.txt"
       if os.path.isfile(list_name):
          os.remove(list_name)

       # use mencoder to generate an animation from the list
       fps = "5" # frames per sec

       # mencoder cannot handle eps files
       if not(variables['eps']):
 
          if (variables['movie']=="mp4"):

             movie_output2 = "Time-Animation.mp4"

             if os.path.isfile(movie_output2):
                os.remove(movie_output2)
             
             string1 = "mencoder mf://@"+list_name+" -o " + movie_output2
             string2 = " -of lavf -lavfopts format=mp4 -ss 1 -ovc "
             string3 = "x264 -x264encopts "
             string4 = "crf=20.0:nocabac:level_idc=30:"
             string5 = "global_header:threads=2 -fps " + fps

             os.system(string1+string2+string3+string4+string5)

          elif (variables['movie']=="avi"):

             movie_output2 = "Time-Animation.avi"

             if os.path.isfile(movie_output2):
                os.remove(movie_output2)

             string1 = "mencoder mf://@"+list_name+" -o " + movie_output2
             string2 = " -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=3000"
             string3 = ":vhq -mf type=png -fps " +fps

             os.system(string1+string2+string3)

          else:
             print "\n---ERROR: unknown movie format---\n"

          os.remove(list_name)

       else:
          print "\n---ERROR: mencoder cannot animate eps images---"
          print   "          list of image files: " + list_name +"\n"


    # get stop time
    stoptime = tm.time()
    duration = (stoptime - starttime)/60.0
    duration = int(1000*duration)/1000.0

    print "--------------------------------------------"

    if (variables['component_plots']):
       print "\nCreated Animation: " + movie_output3

    if (variables['position_plots']):
       print "\nCreated Animation: " + movie_output2

    if (variables['animation']):
       print "\nCreated Animation: " + movie_output

    print "\nImages are located in: "+outputDir +"/"
    print "\nProgram Duration: " + str(duration) + " min"
    print "\n---Complete---\n"


#-----------------------------------------------------------------------------
# Routine for plotting
#
def plotting(variables, x_jump_indices, n, nstep, time, coords):

    # explicitly define colors because we are graphing pieces
    # of the path, and pylab gives each piece a different color
    # by default

    colors = ["r","g","b","y","k","c","m"]
    #col = colors[random.randrange(0,len(colors)-1,1)]
    col = "r"
    if (n+1) % 7 == 1:
        col = colors[0]
    elif (n+1) % 7 == 2:
        col = colors[1]
    elif (n+1) % 7 == 3:
        col = colors[2]
    elif (n+1) % 7 == 4:
        col = colors[3]
    elif (n+1) % 7 == 5:
        col = colors[4]
    elif (n+1) % 7 == 6:
        col = colors[5]
    elif (n+1) % 7 == 0:
        col = colors[6]


    # PLOT X vs T
    # these plots are meant for an animation in time
    if (variables['x_vs_t'] and variables['animation']):

       # unwrap indices for particle n
       # these indices are the places where a particle jumps from 
       # one boundary to another
       ind2 = x_jump_indices[n]
    
       # plot the data, leaving out the jump between boundaries
       if len(ind2) != 0:
          ind = []

          # for each index:
          for num in ind2:

              # only include indices that correspond to current time
              if num <= nstep:
                 ind.append(num)

          # no jumps to worry about (yet)
          if len(ind)==0:
             pylab.plot(time[0:nstep],coords[0,0:nstep],color=col)

          # plot the data in pieces: splitting it where it jumps the
          # boundary, and it jumps the boundary at precisely the
          # the indicies in the ind list
          else:
             pylab.plot(time[0:ind[0]],coords[0,0:ind[0]],color=col)

             for i in range(0,len(ind)-1,1):
                pylab.plot(time[(ind[i]+1):ind[i+1]],
                       coords[0,(ind[i]+1):ind[i+1]], color = col)

             pylab.plot(time[(ind[-1]+1):nstep],
                    coords[0,(ind[-1]+1):nstep],color=col)

       # there are no jumps between boundaries
       else:
          pylab.plot(time[0:nstep],coords[0,0:nstep],color=col)

    # PLOT Y vs T
    # these plots are meant for an animation in time
    elif (variables['y_vs_t'] and variables['animation']):

          pylab.plot(time[0:nstep],coords[1,0:nstep],color=col)

    # PLOT Y vs X
    # these plots are meant for an animation in time
    elif (variables['y_vs_x'] and variables['animation']):

       # unwrap indices for particle n
       ind2 = x_jump_indices[n]
       # plot the ranges, leaving out the jump between boundaries
       if len(ind2) != 0:
          ind = []
          for num in ind2:
              if num <= nstep:
                 ind.append(num)
          if len(ind)==0:
             pylab.plot(coords[0,0:nstep],coords[1,0:nstep],color=col)

          else:
             pylab.plot(coords[0,0:ind[0]],coords[1,0:ind[0]],color=col)

             for i in range(0,len(ind)-1,1):
                pylab.plot(coords[0,(ind[i]+1):ind[i+1]],
                           coords[1,(ind[i]+1):ind[i+1]], color=col)

             pylab.plot(coords[0,(ind[-1]+1):nstep],
                        coords[1,(ind[-1]+1):nstep],color=col)

       else:
          pylab.plot(coords[0,0:nstep],coords[1,0:nstep],color=col)

    # PLOT X vs T position plots
    # plot entire trajectory for each particle
    elif (variables['x_vs_t'] and variables['position_plots']):

       # unwrap indices for particle n
       # the indices correspond to where the particle jumps from 
       # one boundary to another
       ind = x_jump_indices[n]
      
       # plot the data, leaving out the jump between boundaries
       if len(ind) != 0:

             pylab.plot(time[0:ind[0]],coords[0,0:ind[0]],color=col)

             for i in range(0,len(ind)-1,1):
                pylab.plot(time[(ind[i]+1):ind[i+1]],
                       coords[0,(ind[i]+1):ind[i+1]], color = col)

             pylab.plot(time[(ind[-1]+1):],
                    coords[0,(ind[-1]+1):],color=col)

       # no jumps to worry about
       else:
          pylab.plot(time[0:],coords[0,0:],color=col)


    # PLOT Y vs T position plots
    # plot entire trajectory for each particle
    elif (variables['y_vs_t'] and variables['position_plots']):

          pylab.plot(time[0:],coords[1,0:],color=col)


    # PLOT Y vs X position plots
    # plot entire trajectory for each particle
    elif (variables['y_vs_x'] and variables['position_plots']):

       # unwrap indices for particle n
       ind = x_jump_indices[n]
       # plot the ranges, leaving out the jump between boundaries
       if len(ind) != 0:

             pylab.plot(coords[0,0:ind[0]],coords[1,0:ind[0]],color=col)

             for i in range(0,len(ind)-1,1):
                pylab.plot(coords[0,(ind[i]+1):ind[i+1]],
                           coords[1,(ind[i]+1):ind[i+1]],color=col)

             pylab.plot(coords[0,(ind[-1]+1):],
                        coords[1,(ind[-1]+1):],color=col)

       else:
          pylab.plot(coords[0,:],coords[1,:],color=col)



#-----------------------------------------------------------------------------
# This allows for an inputs file to define the variables
#
def read_from_input_file(Variables, Filename):

    if (os.path.isfile(Filename)):

        f = open(Filename)

        for line in f:

            # ignore blanks and lines that start with "#"
            if (not(line.strip().startswith("#") or line.isspace())):

                # get rid of spaces within a line
                line=line.replace(" ","")

                # split the line at "="
                lin = line.split("=",1)

                # loop over all input list variables that are defined 
                for key, value in Variables.iteritems():

                    # variable appears in the input file
                    if (key==lin[0].strip()):

                        # value = variant of true, set to logical true
                        if ((lin[1].strip() == "True") or 
                            (lin[1].strip() == "true") or
                            (lin[1].strip() == "T") or
                            (lin[1].strip() == "t")):
                           Variables[key] = True
                        # value = variant of false, set to logical false
                        elif ((lin[1].strip() == "False") or
                              (lin[1].strip() == "false") or
                              (lin[1].strip() == "f") or
                              (lin[1].strip() == "F")):
                           Variables[key] = False

                        # set variable to value in input file
                        else:
                           Variables[key] = lin[1].strip()

        f.close()

    # input file does not exist
    else:
        print "\n---ERROR: Input file not found---\n"
        sys.exit(2)


#-----------------------------------------------------------------------------
if __name__== "__main__":


    if (len(sys.argv) == 1):
        print "\n---ERROR: No particle data files specified---\n"
        print   "Usage:\n"
        print   "   ./plotparticles.py particle-file1 [particle-file2 ...]\n"
        sys.exit(2)


    main(sys.argv[1:])
        
