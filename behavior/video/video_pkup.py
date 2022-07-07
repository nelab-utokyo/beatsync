import cv2
import os
import pandas as pd

def rewrite_video(datapath, videoname, videoname2, st, en):
    size = (640, 480)
    #fps = 20.1545
    fps = 20.0349
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    video = cv2.VideoWriter(datapath + videoname2, fourcc, fps, size)
    cap = cv2.VideoCapture(datapath + videoname)

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


datapath = "/Users/yoshiki/research/Tlab/behavior/acceleration/20201130/"
dirs = os.listdir(datapath)
dirs = [i for i in dirs if i.startswith("video")]
vdirs = sorted(dirs)
print(dirs)

dirs = sorted(os.listdir(datapath))
dirs = [dirs[i] for i in range(len(dirs)) if dirs[i].startswith("vtime")]
csvs = []
for j in vdirs:
    vfile = j.split(".")[0].split("_",1)[1]
    csvs.extend([i for i in dirs if vfile in i ])


print(csvs)

chose = 0
#rat 6
st = 14161
en = 15587

#rat 8
#st = 6395 - 100
#en = 7617 + 100

#rat 4
#st = 7275
#en = 8735

df = pd.read_csv(datapath + csvs[chose])
print("csv length",len(df["framenum"]))

vname2 = vdirs[chose].rsplit(".",1)[0] + "_pkup.avi"
rewrite_video(datapath, vdirs[chose],vname2, st, en)


