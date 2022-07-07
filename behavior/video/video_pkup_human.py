import cv2
import os
import pandas as pd
import numpy as np
from datetime import datetime as dt

def rewrite_video(datapath, videoname, videoname2, st, en, fps):
    cap = cv2.VideoCapture(datapath + videoname)
    w, h = cap.get(cv2.CAP_PROP_FRAME_WIDTH), cap.get(cv2.CAP_PROP_FRAME_HEIGHT)
    size = (int(w),int(h))
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(datapath + videoname2, fourcc, fps, size)

    if not cap.isOpened():
        return

    counter = 0
    while True:
        counter += 1
        ret, frame = cap.read()
        if counter >= st and counter <= en:
            video.write(frame)
            if counter == en:
                return

patha = "/Users/TP/research/Tlab/behavior/human/20220502/"
pathb = "/Users/TP/research/Tlab/behavior/human/20220513/"

dr = os.listdir(patha)
vdirs = sorted([d for d in dr if d.endswith("mp4") and "pkup" not in d])
drb = os.listdir(pathb)
drb2 = sorted([d for d in drb if d.endswith("mp4") and "pkup" not in d])

vdirs.extend(drb2)

dr = os.listdir(patha)
dirs = sorted([d for d in dr if d.startswith("vtime")])
drb = os.listdir(pathb)
drb2 = sorted([d for d in drb if d.startswith("vtime")])

dirs.extend(drb2)

csvs = []
for j in vdirs:
    vfile = j.split(".")[0].split("_",1)[1]
    csvs.extend([i for i in dirs if vfile in i ])

print(csvs)


chose = 6

datapath = patha if chose < 3 else pathb

df = pd.read_csv(datapath + csvs[chose])
tmp = np.where(df["kinds"] == 1)[0]
pos = np.where(np.diff(tmp) > 5)[0]

st = tmp[0]
en = tmp[pos][0]

print(st, en, len(pos))
t1 = dt.strptime(df["time"][st], '%Y-%m-%d %H:%M:%S.%f')
t2 = dt.strptime(df["time"][en], '%Y-%m-%d %H:%M:%S.%f')
dur = t2 - t1
fps = (en - st)/dur.total_seconds()
print(fps)
print("csv length",len(df["framenum"]))

vname2 = vdirs[chose].rsplit(".",1)[0] + "_pkup.mp4"
rewrite_video(datapath, vdirs[chose],vname2, st, en, fps)


