import subprocess
from os import walk
DIR = "/media/dnx/DnX Data/Dropbox/Education/USP/master2/Grafica/pbrt-importance-sampling/pbrt-v2-master/outputs"
FOLDERS = ["01DP"]
#, "02VW", "03RD", "04PD", "05RW", "06OP", "07LA"
for folder in FOLDERS:
    f = []
    path = DIR+"/"+folder+"/pbrt/"
    for (dirpath, dirnames, filenames) in walk(path):
        f.extend(filenames)
        break
    for fname in f:
        output = subprocess.check_output(['../src/bin/pbrt', path+fname])
        print path+fname
