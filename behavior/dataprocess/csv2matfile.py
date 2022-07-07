import pandas as pd
import numpy as np
import os
import datetime
import scipy.io

def loadfiles(path, dirs, pkup = 0):
    filedir = dirs[pkup]
    flag = False
    fpath = path + filedir + "/"
    files = [d for d in os.listdir(fpath) if d.startswith("AppTag")]
    files = sorted(files)
    if "withlabel" in os.listdir(path):
        vfiles = [d for d in os.listdir(path+"withlabel") if d.startswith("vtime") and filedir.rsplit("_",1)[1] in d and "edited" in d]
        if len(vfiles) > 0:
            flag = True
            df = pd.read_csv(path +"withlabel/"+ vfiles[0], index_col=0)
            print(vfiles)
            df['label'] = pd.to_numeric(df['label'], errors='coerce')
            ls = df['label'].tolist()
            labels = np.array([0 if np.isnan(i) else int(i) for i in ls])
    vfiles = [d for d in os.listdir(path) if d.startswith("vtime") and filedir.rsplit("_",1)[1] in d][0]
    print(vfiles)
    stamp = []
    for j in range(len(files)):
        df = pd.read_csv(fpath + files[j], index_col=0)
        ax_ = df["AccelerationX"]
        ay_ = df["AccelerationY"]
        az_ = df["AccelerationZ"]
        temp = np.array([np.float64(ax_),np.float64(ay_),np.float64(az_)])
        stamp_ = [datetime.datetime.strptime(i, '\t%Y/%m/%d %H:%M:%S.%f') for i in df.index if len(str(i))>4]
        if j == 0:
            acc = np.array(temp)
            time0 = stamp_[0]
        else:
            acc = np.hstack((acc, temp))
        stamp_ = [(st - time0).total_seconds() for st in stamp_]
        xs = np.arange(0,len(stamp_)*10,10)
        xs2 = np.arange(len(df.index))
        stamp.extend(np.interp(xs2,xs,stamp_))
    df = pd.read_csv(path + vfiles, index_col=0)
    stamp = np.array(stamp)
    vstamp = [datetime.datetime.strptime(i, '%Y-%m-%d %H:%M:%S.%f') for i in df["time"] if len(str(i))>4]
    vstamp = [(st - time0).total_seconds() for st in vstamp]
    monoff = np.int8(df["music"] != -1)
    monoff_ = df["music"]
    mkinds = np.int8(df["kinds"])
    ontime = [(vstamp[i+2],monoff_[i+2], mkinds[i+2]) for i in range(len(monoff[:-2])) if monoff[i+2] - monoff[i+1] == 1 and monoff[i] == 0]
    offtime = [(vstamp[i],monoff_[i],mkinds[i]) for i in range(len(monoff[:-2])) if monoff[i+1] - monoff[i] == -1 and monoff[i+2] == 0]
    offtime.append((np.array(vstamp)[-1],np.array(monoff_)[-1],np.array(mkinds)[-1]))
    print(len(ontime), len(offtime))
    click = np.where(df["click"] != 0)[0]
    click = np.vstack((np.array(vstamp)[click],np.array(mkinds)[click],np.array(monoff_)[click])).T
    if flag:
        return acc, stamp, vstamp, ontime, offtime, click, labels
    else:
        return acc, stamp, vstamp, ontime, offtime, click



path = "../../../../acceleration/zemi/20201118/"

dirs = [d for d in os.listdir(path) if d.startswith("AppTag") and os.path.isdir(path+d)]
dirs = sorted(dirs)
print(dirs)

names = ["acc", "stamp", "vstamp", "ontime", "offtime", "click", "labels"]

for i in range(len(dirs)):
    data_ = loadfiles(path, dirs, pkup = i)
    data = {}
    for j in range(len(data_)):
        data[names[j]] = data_[j]
    scipy.io.savemat(path + dirs[i] + ".mat", data)



