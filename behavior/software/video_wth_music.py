import numpy as np
import cv2
import time,datetime
import sys
import os
import pyaudio
import wave
import pandas as pd
import subprocess
import pygame
import pickle

import mmap


def record(frame,datime, framenum, onoff, kinds, click):
    global datafile
    video.write(frame)
    click2 = -int(click == 1)
    datafile.append([framenum, datime, onoff, kinds, click2])

def fin():
    print("finish")
    df = pd.DataFrame(np.array(datafile),columns=column1)
    df.to_csv(filename)
    wf = wave.open(soundname, 'wb')
    wf.setnchannels(CHANNELS)
    wf.setsampwidth(p.get_sample_size(FORMAT))  #width=2 ; 16bit
    wf.setframerate(RATE)
    wf.writeframes(b''.join(sdata))
    wf.close()
    cmd = "ffmpeg -i {} -i {} -c:v copy -c:a copy -b:v 6144k -b:a 448k {}".format(soundname,videoname,videoname2)
    subprocess.call(cmd, shell=True)

def music_play(i,m):
    pygame.mixer.music.load(soundpath + musiclist[m])
    pygame.mixer.music.play(i-1)


#music_setting
musiclist = ["sonata_s75.wav","Mozart_k448.wav","sonata_s200.wav","sonata_s400.wav"]
abbname = "4spNoclick"
rpnum = [1,1,2,4]
soundpath = "C:/Users/Behavior/Desktop/behavior/accelerator/code/sound/"
pygame.mixer.init(frequency=96000, size=-16, channels=2, buffer=4096)


#video recording setting
datetime_now = datetime.datetime.now()
print(datetime_now)
datapath = "C:/Users/Behavior/Desktop/behavior/accelerator/data/{}{:02}{:02}/".format(datetime_now.year,datetime_now.month,datetime_now.day)
os.makedirs(datapath, exist_ok=True)
videoname = datapath+"temp.avi"
videoname2 = datapath+"video_{}_{}{:02}{:02}_{:02}-{:02}.avi".format(abbname,datetime_now.year,datetime_now.month,datetime_now.day,datetime_now.hour,datetime_now.minute)
size = (640, 480)
fps = 24
fourcc = cv2.VideoWriter_fourcc(*'XVID')
video = cv2.VideoWriter(videoname, fourcc, fps, size)

#sound recording setting
FORMAT = pyaudio.paInt16
CHANNELS = 1  #monoral
RATE = 24000
CHUNK = 1024
p = pyaudio.PyAudio()
stream = p.open(format=FORMAT,channels=CHANNELS,rate=RATE,input=True,frames_per_buffer=CHUNK,output = True)
sdata = []
soundname = datapath+"temp.wav"

#csv Output
datafile = []
column1 = ["framenum","time","music","kinds","click"]
filename = datapath+"vtime_{}_{}{:02}{:02}_{:02}-{:02}.csv".format(abbname,datetime_now.year,datetime_now.month,datetime_now.day,datetime_now.hour,datetime_now.minute)

videocapture_num = 0
capture = cv2.VideoCapture(videocapture_num)
cv2.namedWindow("frame")

mode_play_or_stop = 0
framenum = -1
framenum2 = 0
rnd = -1
mnum = np.zeros(len(musiclist))
monoff = 1

last = 8    #the number of music session
offlen = 60 #[s]

t0 = time.time()
t1 = time.time()

while True:
    ret, frame = capture.read()
    datime = str(datetime.datetime.now())
    cv2.imshow("frame", frame)
    c = cv2.waitKey(1)
    if c == 27: # Escキー
        break
    elif c == 0x73:# sキー再生・停止の制御
        if mode_play_or_stop == 0:#停止中の場合　音楽再生モードに切り替え
            mode_play_or_stop = 1
            t1 = time.time()
            print("start", flush=True)
    t = time.time()
    pos = pygame.mixer.music.get_pos()
    if int(pos) == -1 and mode_play_or_stop == 1:
        if monoff == 0 and int(t-t1) >= offlen:
            framenum2 += 1
            monoff = 1
            flag = True
            tmp = np.where(mnum == np.min(mnum))[0]
            rnd = int(np.random.choice(tmp,1))
            print("on", int(framenum2), musiclist[rnd], flush=True)
            mnum[rnd] += 1
            music_play(int(rpnum[rnd]), rnd)
        elif monoff == 1:
            if rnd == 2 and flag:
                flag = False
                music_play(int(rpnum[rnd]), rnd)
            else:
                monoff = 0
                rnd = -1
                t1 = time.time()
                print("off", int(framenum2)+1, flush=True)
                if framenum2 >= last:
                    break
    input = stream.read(CHUNK)
    sdata.append(input)
    if framenum != int((t-t0)*fps):
        framenum = int((t-t0)*fps)
        conoff = 0
        record(frame, datime, framenum, pos, rnd, conoff)


capture.release()
cv2.destroyAllWindows()
video.release()

fin()
