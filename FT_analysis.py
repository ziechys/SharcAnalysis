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
import os, sys, shutil, re, datetime, readline
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
except:
    print('Plotting not possible (probably because pylab/matplotlib is not installed or because there is no X connection)')
    plot_possible = False

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
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.trajana_nma')

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
    print('')
  else:
    INFOS['plot']=False

  print('')
  intervallist=[]
  print(centerstring('Time steps to be analysed',60,'-'))
  print('\nPlease enter the time step intervals for which the statistical analysis should be carried out. ')
  print('')
  while True:
    interval=question('Time step interval: ',int,[0,maxlen])
    #endtime=question('End time of interval: ',int,[maxlen])[0]
    print('')
    #interval=[starttime,endtime]
    intervallist.append(interval)
    moreinterval=question('Do you want to add another time interval for analysis?',bool,False)
    if not moreinterval:
      break
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
  print(centerstring('Results directory',60,'-'))
  print('Please give the name of the subdirectory to be used for the results (use to save similar analysis in separate subdirectories).')
  destin=question('Name for subdirectory?',str,'nma')
  INFOS['descr']=destin
  print('')

  return INFOS

# ======================================================================= # 
def ft_analysis(INFOS):
    """
    FT analysis analysis. Typically this script is carried out.
    """
    print('Preparing FT analysis ...')

    num_steps=INFOS['numsteps']
    descr=INFOS['descr']
    out_dir = os.path.join('NMA',descr)
    ref_struc_file=INFOS['refstruc']
    ref_struc_type=INFOS['reftype']
    ana_ints=INFOS['interval']
    vibration_file=INFOS['refvib']
    dt=INFOS['timestep']
    abs_list=INFOS['symmmodes']
    neg_list=INFOS['negmodes']
    plot=INFOS['plot']
    mawe=INFOS['massweight']    

    try:
        os.makedirs(out_dir)
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

    for i in range(ntraj):  #loop over all trajectories

      print('Reading trajectory ' + str(files[i]) + ' ...')
      folder_name = str(files[i])[:-10]

      ################################################
      #
      # Analysis A: Create geometry file based on SHARC's geo.py and do individual FT analysis

      os.system("python3 geo.py -g " + folder_name + "/output.xyz -t " + INFOS['timestep'] + " < " + INFOS['geo'] + " > " + folder_name + "/FT_Geo.out")
      sys.exit(0)

      ################################################
      #
      # Analysis B: Create coherent structure by adding up all structures at all time steps (see essential dynamics)

      #call traj_manip class for single traj: Creates a trajectory object out of Nt structure objects
      trajectory = traj_manip.trajectory(folder_name, ref_struc, dt=dt,GSstop=INFOS['GShop'])

      
    
    
    print('Data processing finished.')
    
    #if plot: plot_summary(INFOS)
    
def plot_summary(INFOS):
    # plotting

    descr=INFOS['descr']
    out_dir = os.path.join('NMA',descr)
    num_steps=INFOS['numsteps']

    if plot_possible:
        print('Drawing plots ...')
        
        # Plots for time dependent cross averages
        plotting.mean_std_from_files(mean_file=out_dir+'/mean_against_time.txt',out_dir=out_dir+'/time_plots',xlist=[INFOS['timestep'] * i for i in range(num_steps)],std_file=out_dir+'/std_against_time.txt')
        
        # Bar graphs with the standard deviation of time dependent cross averages
        plotting.bars_from_file(in_file=out_dir+'/total_av.txt', out_dir=out_dir+'/bar_graphs/total_av')
        plotting.bars_from_file(in_file=out_dir+'/total_std.txt', out_dir=out_dir+'/bar_graphs/total_std')
        plotting.bars_from_file(in_file=out_dir+'/total_weighted_av.txt', out_dir=out_dir+'/bar_graphs/total_weighted_av')
        plotting.bars_from_file(in_file=out_dir+'/total_weighted_std.txt', out_dir=out_dir+'/bar_graphs/total_weighted_std')
        plotting.bars_from_file(in_file=out_dir+'/coh_av.txt', out_dir=out_dir+'/bar_graphs/coh_av')
        plotting.bars_from_file(in_file=out_dir+'/coh_std.txt', out_dir=out_dir+'/bar_graphs/coh_std')
        plotting.bars_from_file(in_file=out_dir+'/coh_weighted_av.txt', out_dir=out_dir+'/bar_graphs/coh_weighted_av')
        plotting.bars_from_file(in_file=out_dir+'/coh_weighted_std.txt', out_dir=out_dir+'/bar_graphs/coh_weighted_std')
    else:
        print('Plotting not possible')
    
def main():
    displaywelcome()
    open_keystrokes()

    INFOS=get_general()
    #ft_analysis(INFOS)
                    
    close_keystrokes()
                
if __name__=='__main__':
  try:
    main()
  except KeyboardInterrupt:
    print('\nCtrl+C makes me a sad SHARC ;-(\n')
    quit(0)

#print modes_ind, plot_list, modes_list
