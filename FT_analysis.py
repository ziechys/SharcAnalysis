#!/usr/bin/env python3

#******************************************
#
#   SHARC Analysis script
#   Written in style with original SHARC scripts and using some original SHARC classes
#   No affiliation with the original SHARC distributor
#   No guarantee for functionality
#
#******************************************

"""
author: Karl Michael Ziems
version: 1.0
descr: Analyse structural dynamics of coherent classical motion during relaxation dynamics by means of Fourier Transformation to extract frequency data
"""

#import warnings
#warnings.filterwarnings("error")
#warnings.simplefilter("ignore", category=ResourceWarning) 
import os, sys, shutil, re, datetime, readline, copy
from shutil import copyfile
os.chdir(os.path.dirname(sys.argv[0]))
sys.path.append('sharclib')
try:
    import numpy
except ImportError:
    print('numpy package not installed')
    sys.exit()
try:
    import file_handler, vib_molden, traj_manip, struc_linalg
except ImportError:
    print('file_handler, vib_molden, traj_manip or struc_linalg not found. They should be part of this package. Check the installation and if $SHARC/../lib is part of the PYTHONPATH environment variable.')
    sys.exit()
try:
    import plotting
    plot_possible = True
    import matplotlib
    import matplotlib.pyplot as plt
except:
    print('Plotting not possible (probably because pylab/matplotlib is not installed or because there is no X connection)')
    plot_possible = False

from scipy.interpolate import CubicSpline
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.fft import fft, ifft, fftfreq, rfft, fftshift
from scipy import optimize
import math

try:
    from openbabel import openbabel
except ImportError:
    print('openbabel.py package not installed')
    sys.exit()
version='1.0'
versiondate=datetime.date(2021,6,4)

# ======================================================================= #
def centerstring(string,n,pad=' '):
  l=len(string)
  if l>=n:
    return string
  else:
    return  pad*int((n-l+1)/2)+string+pad*int((n-l)/2)

# ======================================================================= #
def displaywelcome():
  print('Script for structural dynamics Fourier tranformation analysis started ...\n')
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('FT analysis for SHARC dynamics',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Karl Michael Ziems',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='''
This script reads output.xyz files, create a coherent motion, extracts geometry data and performs FT analysis of coherent motion and single trajectories.
  '''
  print(string)

# ======================================================================= #
def readfile(filename):
  try:
    f=open(filename)
    out=f.readlines()
    f.close()
  except IOError:
    print('File %s does not exist!' % (filename))
    sys.exit(12)
  return out

# ======================================================================= #
def writefile(filename,content):
  # content can be either a string or a list of strings
  try:
    f=open(filename,'w')
    if isinstance(content,list):
      for line in content:
        f.write(line)
    elif isinstance(content,str):
      f.write(content)
    else:
      print('Content %s cannot be written to file!' % (content))
    f.close()
  except IOError:
    print('Could not write to file %s!' % (filename))
    sys.exit(13)
# ======================================================================= #

def open_keystrokes():
  global KEYSTROKES
  KEYSTROKES=open('KEYSTROKES.tmp','w')

# ======================================================================= #
def close_keystrokes():
  KEYSTROKES.close()
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.ft')

# ======================================================================= #
def question(question,typefunc,default=None,autocomplete=True,ranges=False):
  if typefunc==int or typefunc==float:
    if not default==None and not isinstance(default,list):
      print('Default to int or float question must be list!')
      quit(1)
  if typefunc==str and autocomplete:
    readline.set_completer_delims(' \t\n;')
    readline.parse_and_bind("tab: complete")    # activate autocomplete
  else:
    readline.parse_and_bind("tab: ")            # deactivate autocomplete

  while True:
    s=question
    if default!=None:
      if typefunc==bool or typefunc==str:
        s+= ' [%s]' % (str(default))
      elif typefunc==int or typefunc==float:
        s+= ' ['
        for i in default:
          s+=str(i)+' '
        s=s[:-1]+']'
    if typefunc==str and autocomplete:
      s+=' (autocomplete enabled)'
    if typefunc==int and ranges:
      s+=' (range comprehension enabled)'
    s+=' '

    line=input(s)
    line=re.sub('#.*$','',line).strip()
    if not typefunc==str:
      line=line.lower()

    if line=='' or line=='\n':
      if default!=None:
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return default
      else:
        continue

    if typefunc==bool:
      posresponse=['y','yes','true', 't', 'ja',  'si','yea','yeah','aye','sure','definitely']
      negresponse=['n','no', 'false', 'f', 'nein', 'nope']
      if line in posresponse:
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return True
      elif line in negresponse:
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return False
      else:
        print('I didn''t understand you.')
        continue

    if typefunc==str:
      KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
      return line

    if typefunc==float:
      # float will be returned as a list
      f=line.split()
      try:
        for i in range(len(f)):
          f[i]=typefunc(f[i])
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return f
      except ValueError:
        print('Please enter floats!')
        continue

    if typefunc==int:
      # int will be returned as a list
      f=line.split()
      out=[]
      try:
        for i in f:
          if ranges and '~' in i:
            q=i.split('~')
            for j in range(int(q[0]),int(q[1])+1):
              out.append(j)
          else:
            out.append(int(i))
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return out
      except ValueError:
        if ranges:
          print('Please enter integers or ranges of integers (e.g. "-3~-1  2  5~7")!')
        else:
          print('Please enter integers!')
        continue

# ======================================================================= #
def get_general():
  """
  Parsing information from terminal
  """

  INFOS={}

  print(centerstring('Paths to trajectories',60,'-'))
  print('\nPlease enter the paths to all directories containing the "TRAJ_0XXXX" directories.\nE.g. Sing_2/ and Sing_3/. \nPlease enter one path at a time, and type "end" to finish the list.')
  count=0
  paths=[]
  while True:
    path=question('Path: ',str,'end')
    if path=='end':
      if len(paths)==0:
        print('No path yet!')
        continue
      print('')
      break
    path=os.path.expanduser(os.path.expandvars(path))
    if not os.path.isdir(path):
      print('Does not exist or is not a directory: %s' % (path))
      continue
    if path in paths:
      print('Already included.')
      continue
    ls=os.listdir(path)
    print(ls)
    for i in ls:
      if 'TRAJ' in i:
        count+=1
    print('Found %i subdirectories in total.\n' % count)
    paths.append(path)
  INFOS['paths']=paths
  print('Total number of subdirectories: %i\n' % (count))

  # try to obtain the maximum number of time steps in the trajectories
  maxlen=0
  dt=0.0
  forbidden=['crashed','running','dead','dont_analyze']
  for idir in INFOS['paths']:
    ls=os.listdir(idir)
    for itraj in ls:
      if not 'TRAJ_' in itraj:
        continue
      path=idir+'/'+itraj
      pathfile=path+'/output.lis'
      if not os.path.isfile(pathfile):
        continue
      lstraj=os.listdir(path)
      valid=True
      for i in lstraj:
        if i.lower() in forbidden:
          valid=False
          break
      if not valid:
        continue
      f=readfile(pathfile)
      for line in f:
        if '#' in line:
          continue
        s=line.split()
        step=int(s[0])
        if step>maxlen:
          maxlen=step
        if dt==0.:
          dt=float(s[1])

  print(centerstring('Path to normal mode file',60,'-'))
  print('\nPlease enter the path to the Molden normal mode file for your molecule. The contained geometry will be used as reference geometry.\n(Atomic order must be the same as in the trajectories!)')
  print('')
  refvib=question('Path: ',str)
  refvib=os.path.expanduser(os.path.expandvars(refvib))
  INFOS['refvib']=refvib
  INFOS['refstruc']=refvib
  INFOS['reftype']='molden'
  print('')
  
  INFOS['massweight']=question('Do you wish to use mass weighted normal modes?',bool,True)
  print('')

  print(centerstring('Number of total steps in your trajectories',60,'-'))
  #print '\n total simulation time *2 and +1 if timestep is 0.5 fs'
  print('')
  while True:
    numsteps=question('Number of time steps: ',int,[maxlen+1])[0]
    if numsteps <=0:
      print('Number of steps must be positive!')
      continue
    break
  INFOS['numsteps']=numsteps
  print('')

  print(centerstring('The time step of your calculation',60,'-'))
  print('')
  while True:
    timestep=question('Length of time step: ',float,[dt])[0]
    if timestep <=0.:
      print('Time step must be positive!')
      continue
    break
  INFOS['timestep']=timestep
  print('')

  if plot_possible:
    print(centerstring('Automatic plot creation',60,'-'))
    print('')
    autoplot=question('Do you want to automatically create plots of your data?',bool,False)
    INFOS['plot']=autoplot
    if autoplot:
      INFOS['labels']=question('Specify labels for legend based on internal coordinates used in FT analysis (comma separated)',str,autocomplete=False).split(',')
      INFOS['limit']=question('Specify x axis range (in cm^-1)',float,[0,4000])
      INFOS['live']=question('Do you want to live analyse via matplotlib graphic interface?',bool,False)
      INFOS['onefolder']=question('Do you want to save all graphs in one directory?',bool,True)
      if INFOS['onefolder']:
        INFOS['savedir']=question('Name for directory?',str,'FT')
  else:
    INFOS['plot']=False

  print('')
  intervallist=[]
  print(centerstring('Time steps to be analysed',60,'-'))
  print('\nPlease enter the time step intervals for which the FT analysis should be carried out. ')
  print('')
#ToDo: Implement several intervals
#  while True:
#    interval=question('Time step interval: ',int,[0,maxlen])
#    #endtime=question('End time of interval: ',int,[maxlen])[0]
#    print('')
#    #interval=[starttime,endtime]
#    intervallist.append(interval)
#    moreinterval=question('Do you want to add another time interval for analysis?',bool,False)
#    if not moreinterval:
#      break
  interval=question('Time step interval: ',int,[0,maxlen])
  intervallist.append(interval)
  INFOS['interval']=intervallist
  print('')
  GShop=question('Do you want to only regard data until GShop. Feature only available for forced GS hop dynamics.',bool,False)
  INFOS['GShop'] = GShop
  print('')
  print(centerstring('Geometry file',60,'-'))
  print('Please give the geometry file indicating the internal coordinates for FT analysis. The file should be given in SHARC geo.py format!')
  geoin=question('Name for geometry file?',str,'Geo.inp')
  INFOS['geo'] = geoin
  print('')
  print(centerstring('Post-processing',60,'-'))
  post=question('Do you want the FT analysis date to be post-processed for improved resolution?',bool,True)
  INFOS['post'] = post
  print('')
  print(centerstring('Detailed Analysis',60,'-'))
  print('This script will create a coherent wave packet out of all trajectory runs and anlyse the structural dynamics via FT')
  INFOS['all']=question('Do you want to additonally analyse each individual trajectory?',bool,True)
  print('')

  return INFOS

def do_FT(time, data,plottime=False,damping=1.0,padding=0.0,mean=False,corrdamp=False):
  """
  FT anaylsis of data. 
  Possible FT options are: splines, damping, padding, mean value substraction, post-damping correction.
  """
  #Some useful conversion
  cm2au = 4.5563352529120*10**(-6)
  au2fs = 0.024188843265857

  #Process data: Mean offset and damping
  if mean:
      offset = numpy.mean(data[:   round(math.floor(len(data)/numpy.pi/2) * numpy.pi * 2 )   ])
      #offset = numpy.mean(data)
      print('Data offset found of ', offset)
      data = data - offset
  data = data * damping

  if plottime:
      #Plot input data in time
      plt.plot(time*au2fs, data,label=legend)
      plt.xlabel('t / fs')
      plt.ylabel('E / eV')
      plt.show()
  
  #FFT to energy/freq
  if padding > 0.0:
      gridpoints = padding
  else:
      gridpoints = len(data)
  dt = abs(time[1] - time[0])
  dE = (2*numpy.pi)/(dt)
  FFT_energy = rfft(data,n=gridpoints)
  energy = numpy.linspace(0, dE/2, num=gridpoints//2+1, endpoint=True)
  #Plot FFT results
  if corrdamp:
      FFT_damp = rfft(damping,n=gridpoints)
      energy_out = energy/cm2au
      data_out   = abs(FFT_energy/FFT_damp)/numpy.max( abs(FFT_energy/FFT_damp))
  else:
      energy_out = energy[1:]/cm2au #remove first point cause heavily influenced by damping
      data_out   = abs(FFT_energy[1:])/numpy.max( abs(FFT_energy[1:]))

  return numpy.array(energy_out),numpy.array(data_out)

# ======================================================================= # 
def ft_analysis(INFOS):
    """
    FT analysis analysis. Typically this script is carried out.
    """

    au2fs = 0.024188843265857

    print('Preparing FT analysis ...')

    num_steps=INFOS['numsteps']
    ref_struc_file=INFOS['refstruc']
    ref_struc_type=INFOS['reftype']
    ana_ints=INFOS['interval']
    vibration_file=INFOS['refvib']
    dt=INFOS['timestep']
    plot=INFOS['plot']
    mawe=INFOS['massweight']    

    if plot:
      if INFOS['onefolder']:
        try:
            os.makedirs(INFOS['savedir'])
        except OSError:
            print('Output directory could not be created. It either already exists or you do not have writing access.')
    
    # define reference structure from molden file
    ref_struc = struc_linalg.structure('ref_struc') 
    ref_struc.read_file(ref_struc_file, ref_struc_type)
    num_at = ref_struc.ret_num_at()
    mol_calc = struc_linalg.mol_calc(def_file_path=ref_struc_file, file_type=ref_struc_type)

    forbidden=['crashed','running','dead','dont_analyze']
    width=30
    files=[]
    ntraj=0
    print('Checking the directories...')
    for idir in INFOS['paths']:
      ls=os.listdir(idir)
      for itraj in ls:
        if not 'TRAJ_' in itraj:
          continue
        path=idir+'/'+itraj
        s=path+' '*(width-len(path))
        pathfile=path+'/output.xyz'
        if not os.path.isfile(pathfile):
          s+='%s NOT FOUND' % (pathfile)
          print(s)
          continue
        lstraj=os.listdir(path)
        valid=True
        for i in lstraj:
          if i.lower() in forbidden:
            s+='DETECTED FILE %s' % (i.lower())
            print(s)
            valid=False
            break
        if not valid:
          continue
        s+='OK'
        print(s)
        ntraj+=1
        files.append(pathfile)
    print('Number of trajectories: %i' % (ntraj))
    if ntraj==0:
      print('No valid trajectories found, exiting...')
      sys.exit(0)

    # define time-resolved arrays
    cross_num_array = numpy.zeros(num_steps)
    cross_sum_array = numpy.zeros([num_steps,num_at*3], float) # a number for every time step and coordinate; sum, has to be divided by num_array?

    first = True

    for i in range(ntraj):  #loop over all trajectories

      print('Reading trajectory ' + str(files[i]) + ' ...')
      folder_name = str(files[i])[:-10]

      ################################################
      #
      # Analysis A: Create geometry file based on SHARC's geo.py and do individual FT analysis

      #call traj_manip class for single traj: Creates a trajectory object out of Nt structure objects
      trajectory = traj_manip.trajectory(folder_name, ref_struc, dt=dt,GSstop=INFOS['GShop'])

      if INFOS['all']:
        print('Performing internal coordinate analysis of individual trajectory...')
        #ToDo: Change to python3 for Linux machine
        os.system("python geo.py -g " + folder_name + "/output.xyz -t " + str(INFOS['timestep']) + " < " + INFOS['geo'] + " > " + folder_name + "/FT_Geo.out")
       
        # Read Geo data from specific internal coordinates
        time, *geomdata = numpy.loadtxt(folder_name+'/FT_Geo.out',unpack=True)
        time = time/au2fs

        #Restrict to excited state dynamics if requested
        if INFOS['GShop']:
          time = time[:trajectory.hop]
          for nr in range(numpy.shape(geomdata)[0]):
            geomdata[nr] = geomdata[nr][:trajectory.hop]

        #Interval restriction
        time = time[INFOS['interval'][0][0]:INFOS['interval'][0][1]]
        for nr in range(numpy.shape(geomdata)[0]):
            geomdata[nr] = geomdata[nr][INFOS['interval'][0][0]:INFOS['interval'][0][1]]

        #Damping function
        damping = numpy.array([numpy.cos( ( i /(time[-1]) ) *  (numpy.pi/2) )**2 for i in time  ])

        print('FT analysis of selected internal degrees of freedom...')

        #Prepare data arrays
        padding_nr = 10000
        xdata = numpy.zeros( (numpy.shape(geomdata)[0],numpy.shape(geomdata)[1]//2) )
        ydata = numpy.zeros( (numpy.shape(geomdata)[0],numpy.shape(geomdata)[1]//2) )
        if INFOS['post']:
          xdata_post = numpy.zeros( (numpy.shape(geomdata)[0],padding_nr//2+1))
          ydata_post = numpy.zeros( (numpy.shape(geomdata)[0],padding_nr//2+1))

        #Do actual FT
        for nr in range(numpy.shape(geomdata)[0]):
          xdata[nr],ydata[nr] = do_FT(time,geomdata[nr])
          if INFOS['post']:
            xdata_post[nr],ydata_post[nr] = do_FT(time,geomdata[nr],mean=True,damping=damping,padding=padding_nr,corrdamp=True)

        #Save data
        header = 'frequency [cm^-1]   '
        if plot:
          for i in INFOS['labels']:
            header += '|      '+ i + '      '
        numpy.savetxt(folder_name + 'FT.out', numpy.transpose([xdata[0],*ydata]),header=header)
        numpy.savetxt(folder_name + 'FT_processed.out', numpy.transpose([xdata_post[0],*ydata_post]),header=header)

        #Add up all FT data
        if INFOS['post']:
          if first:
            #summed_FTdata      = numpy.zeros_like(ydata) Not possible due to different length for raw data - splining ToDo
            summed_FTdata_post = numpy.zeros_like(ydata_post)
            first = False
          #summed_FTdata      = summed_FTdata + ydata
          summed_FTdata_post = summed_FTdata_post + ydata_post

        #Plotting
        if INFOS['plot']:
          for nr in range(numpy.shape(geomdata)[0]):
            plt.plot(xdata[nr],ydata[nr],label=INFOS['labels'][nr])
          plt.xlabel(r'$\omega / cm^{-1}$')
          plt.xlim(INFOS['limit'])
          plt.legend()
          plt.savefig(folder_name + 'FT.pdf')
          if INFOS['live']:
            plt.show()
          plt.close()

          if INFOS['post']:
            for nr in range(numpy.shape(geomdata)[0]):
              plt.plot(xdata_post[nr],ydata_post[nr],label=INFOS['labels'][nr])
            plt.xlabel(r'$\omega / cm^{-1}$')
            plt.xlim(INFOS['limit'])
            plt.legend()
            plt.savefig(folder_name + 'FT_processed.pdf')
            if INFOS['live']:
              plt.show()
            plt.close()

          if INFOS['onefolder']:
            copyfile(folder_name + 'FT.pdf ',INFOS['savedir'] + '/' + folder_name[-11:-1]+'_FT.pdf')
            if INFOS['post']:
              copyfile(folder_name + 'FT_processed.pdf ',INFOS['savedir'] + '/' + folder_name[-11:-1]+'_FT_processed.pdf')
      
      #sys.exit(0)

      ################################################
      #
      # Analysis B: Create coherent structure by adding up all structures at all time steps (see essential dynamics)

      #Returns a 3N x T matrix with all the coordinates of the timesteps (could be done with ref structure - think about it)    
      coor_matrix = trajectory.ret_coor_matrix()
      
      # Obtain time-resovled coherent structure
      for nr, tstep in enumerate(coor_matrix):
        cross_num_array[nr] += 1
        cross_sum_array[nr] += tstep

      #sys.exit(0)
    
    print('')
    print('Coherent time-resolved analysis: Processing data ...')

    #Truncate data to limits with actual trajectories in it (no divsion by zero)
    for ind,num in enumerate(cross_num_array):
      if num == 0:   
        cross_num_array = cross_num_array[0:ind]
        cross_sum_array = cross_sum_array[0:ind]
        num_steps = ind # num_steps has to be passed as an argument. so it can be changed here.
        break

    #Create average time-resolved structure
    cross_mean_array = numpy.zeros_like(cross_sum_array)
    for i in range(len(cross_num_array)):
      cross_mean_array[i] = cross_sum_array[i] / cross_num_array[i]
    #Problem: overweighted by late trajectories
    #Mabye define a maximum number of trajectories that should still contribute. As percentage of all trajs?
    #plt.plot(numpy.linspace(0,num_steps*0.5,num_steps),cross_num_array)
    #plt.show()

    #Reshape to coor matrix format
    cross_mean_array = cross_mean_array.reshape(num_steps,num_at,3)

    #Data files needed
    coh_outdir ='coherent_traj.xyz'
    coh_geo_outdir = 'coh_Geo.out'
    coh_ft_outdir = 'coh_FT.out'
    coh_ft_processed_outdir = 'coh_FT_processed.out'
    coh_pdf = 'coh_FT.pdf'
    coh_processed_pdf = 'coh_FT_processed.pdf'
    incoh_ft_outdir = 'incoh_FT.out'
    incoh_ft_processed_outdir = 'incoh_FT_processed.out'
    incoh_pdf = 'incoh_FT.pdf'
    incoh_processed_pdf = 'incoh_FT_processed.df'
    if plot:
      if INFOS['onefolder']:
        coh_outdir = INFOS['savedir'] + '/coherent_traj.xyz'
        coh_geo_outdir = INFOS['savedir'] + '/coh_FT_Geo.out'
        coh_ft_outdir = INFOS['savedir'] + '/coh_FT.out'
        coh_ft_processed_outdir = INFOS['savedir'] + '/coh_FT_processed.out'
        coh_pdf = INFOS['savedir'] + '/coh_FT.pdf'
        coh_processed_pdf = INFOS['savedir'] + '/coh_FT_processed.pdf'
        incoh_ft_outdir = INFOS['savedir'] + '/incoh_FT.out'
        incoh_ft_processed_outdir = INFOS['savedir'] + '/incoh_FT_processed.out'
        incoh_pdf = INFOS['savedir'] + '/incoh_FT.pdf'
        incoh_processed_pdf = INFOS['savedir'] + '/incoh_FT_processed.pdf'

    #Converting structure in coor matrix format to xyz output
    ref_struc.make_trajectory(cross_mean_array,coh_outdir)

    #Convert to internal coordinates data
    os.system("python geo.py -g " + coh_outdir + " -t " + str(INFOS['timestep']) + " < " + INFOS['geo'] + " > " + coh_geo_outdir) 
    # Read Geo data from specific internal coordinates
    time, *geomdata = numpy.loadtxt(coh_geo_outdir,unpack=True)
    time = time/au2fs

    #Interval restriction
    if INFOS['interval'][0][1] > num_steps:
      print('Interval ending was larger then simulation time.')
      INFOS['interval'][0][1] = num_steps
      print('Interval set to ', INFOS['interval'][0])
    time = time[INFOS['interval'][0][0]:INFOS['interval'][0][1]]
    for nr in range(numpy.shape(geomdata)[0]):
        geomdata[nr] = geomdata[nr][INFOS['interval'][0][0]:INFOS['interval'][0][1]]

    #Damping function
    damping = numpy.array([numpy.cos( ( i /(time[-1]) ) *  (numpy.pi/2) )**2 for i in time  ])

    print('FT analysis of selected internal degrees of freedom...')

    #Prepare data arrays
    padding_nr = 10000
    xdata = numpy.zeros( (numpy.shape(geomdata)[0],numpy.shape(geomdata)[1]//2) )
    ydata = numpy.zeros( (numpy.shape(geomdata)[0],numpy.shape(geomdata)[1]//2) )
    if INFOS['post']:
      xdata_post = numpy.zeros( (numpy.shape(geomdata)[0],padding_nr//2+1))
      ydata_post = numpy.zeros( (numpy.shape(geomdata)[0],padding_nr//2+1))

    #Do actual FT
    for nr in range(numpy.shape(geomdata)[0]):
      xdata[nr],ydata[nr] = do_FT(time,geomdata[nr])
      if INFOS['post']:
        xdata_post[nr],ydata_post[nr] = do_FT(time,geomdata[nr],mean=True,damping=damping,padding=padding_nr,corrdamp=True)

    #Save data
    header = 'frequency [cm^-1]   '
    if plot:
      for i in INFOS['labels']:
        header += '|      '+ i + '      '
    numpy.savetxt( coh_ft_outdir, numpy.transpose([xdata[0],*ydata]),header=header)
    if INFOS['post']:
      numpy.savetxt( coh_ft_processed_outdir, numpy.transpose([xdata_post[0],*ydata_post]),header=header)

    #Plotting
    if INFOS['plot']:
      for nr in range(numpy.shape(geomdata)[0]):
        plt.plot(xdata[nr],ydata[nr],label=INFOS['labels'][nr])
      plt.xlabel(r'$\omega / cm^{-1}$')
      plt.xlim(INFOS['limit'])
      plt.legend()
      plt.savefig(coh_pdf)
      if INFOS['live']:
        plt.show()
      plt.close()

      if INFOS['post']:
        for nr in range(numpy.shape(geomdata)[0]):
          plt.plot(xdata_post[nr],ydata_post[nr],label=INFOS['labels'][nr])
        plt.xlabel(r'$\omega / cm^{-1}$')
        plt.xlim(INFOS['limit'])
        plt.legend()
        plt.savefig(coh_processed_pdf)
        if INFOS['live']:
          plt.show()
        plt.close()


    ################################################
    #
    # Analysis C: Add up all FT results and renormalize

    if INFOS['all'] and INFOS['post']:
      print('')
      print('Create summed-up FT analysis over all trajectores')

      #Normalise to 1
      for i in range(len(summed_FTdata_post)):
        summed_FTdata_post[i] = summed_FTdata_post[i] / numpy.max(summed_FTdata_post[i])

      #Save data
      header = 'frequency [cm^-1]   '
      if plot:
        for i in INFOS['labels']:
          header += '|      '+ i + '      '
      numpy.savetxt( incoh_ft_processed_outdir, numpy.transpose([xdata_post[0],*summed_FTdata_post]),header=header)

      #Plotting
      for nr in range(numpy.shape(geomdata)[0]):
        plt.plot(xdata_post[0],summed_FTdata_post[nr],label=INFOS['labels'][nr])
      plt.xlabel(r'$\omega / cm^{-1}$')
      plt.xlim(INFOS['limit'])
      plt.legend()
      plt.savefig(incoh_processed_pdf)
      if INFOS['live']:
        plt.show()
      plt.close()

    print('Data processing finished.')
    
def main():
    displaywelcome()
    open_keystrokes()

    INFOS=get_general()
    ft_analysis(INFOS)
                    
    close_keystrokes()
                
if __name__=='__main__':
  try:
    main()
  except KeyboardInterrupt:
    print('\nCtrl+C makes me a sad SHARC ;-(\n')
    quit(0)

#print modes_ind, plot_list, modes_list
