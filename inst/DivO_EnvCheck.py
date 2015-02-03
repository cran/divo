#divo EnvCheck module
#Version: 0.1.1
#Autors: Maciej Pietrzak, Michal Seweryn, Grzegorz Rempala
#Maintainer: Maciej Pietrzak <pietrzak.20@osu.edu>
#License: GPL (>=2)

import numpy as np

def nptest():
    test=np.load('WriteTest')
    np.save("divoEnvCheck_out", test+2)

try: nptest()
except: pass

