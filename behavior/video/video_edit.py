import cv2
import os
import pandas as pd

def rewrite_video(datapath, videoname, videoname2, texts):
    size = (640, 480)
    fps = 24
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    video = cv2.VideoWriter(datapath + videoname2, fourcc, fps, size)
    cap = cv2.VideoCapture(datapath + videoname)

    if not cap.isOpened():
        return

    n = 0
    font = cv2.FONT_HERSHEY_SIMPLEX 
    color = (200,200,200)
    org = (5,100)
    font_scale=0.7
    thickness=2
    while True:
        ret, frame = cap.read()
        if ret:
            if n < len(texts):
                cv2.putText(frame, str(texts[n]), org, font, font_scale, color, thickness, cv2.LINE_AA)
            video.write(frame)
            n += 1
        else:
            print("video length",n)
            return

datapath = "/Users/yoshiki/research/Tlab/behavior/acceleration/20201030/"
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

df = pd.read_csv(datapath + csvs[chose])
print("csv length",len(df["framenum"]))

vname2 = vdirs[chose].rsplit(".",1)[0] + "_edited.avi"
rewrite_video(datapath, vdirs[chose],vname2, df["framenum"])


