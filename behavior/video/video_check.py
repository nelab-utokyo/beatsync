import cv2
import os
import pandas as pd

#checking the number of frames
def rewrite_video(datapath, videoname, videoname2): 
    size = (640, 480)
    #fps = 20.1545
    fps = 20.0349
    cap = cv2.VideoCapture(datapath + videoname)

    if not cap.isOpened():
        return

    counter = 0
    while True:
        counter += 1
        ret, frame = cap.read()
        if ret:
            counter += 1
        else:
            print(counter)
            return

#exporting videos with frame number
def rewrite_video2(datapath, videoname, videoname2, texts):
    cap = cv2.VideoCapture(datapath + videoname)
    size = (int(cap.get(3)), int(cap.get(4)))
    print(size)
    fps = 20
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    video = cv2.VideoWriter(datapath + videoname2, fourcc, fps, size)
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


datapath = "/Users/yoshiki/research/Tlab/behavior/acceleration/20201207/"
dirs = os.listdir(datapath)
dirs = [i for i in dirs if i.endswith("pkup.mov")]
vdirs = sorted(dirs)
print(dirs)

vname2 = "1st_pkup_edited.avi"
#rewrite_video(datapath, dirs[0],vname2)
rewrite_video2(datapath, dirs[0],vname2, range(500))



