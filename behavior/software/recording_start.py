import subprocess
from subprocess import Popen
import os
import sys


path = "C:/Users/Behavior/Desktop/behavior/accelerator/code/"
cmd1 = "python {}video_wth_music.py".format(path)

proc = Popen( cmd1,shell=True )

os.chdir(path)
cmd2 = "python Graph.py -t COM4 -a".format(path)

proc2 = Popen( cmd2,shell=True )
