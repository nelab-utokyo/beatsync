import cv2
import os, copy
import numpy as np
import scipy.io


def detect_object(datapath, videoname, videoname2, param):
    size = (param[3]-param[2],param[1]-param[0])
    fps = 30
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(datapath + videoname2 + ".mp4", fourcc, fps, size)
    cap = cv2.VideoCapture(datapath + videoname)
    pos = []
    count = 0
    while True:
        ret, frame = cap.read()
        if ret:
            frame_ = frame[param[0]:param[1],param[2]:param[3],:]
            tmpframe_ = np.array(frame_, dtype = "uint8")
            tmpframe_[:,:,2] = np.where(tmpframe_[:,:,2] < param[4], 0, tmpframe_[:,:,2]) #video2
            tmpframe_[:,:,2] = np.uint8(tmpframe_[:,:,2]/param[5]) #video1

            tmpframe_ = np.argmax(tmpframe_, axis=2)
            tmpframe = np.uint8((tmpframe_ == 2)*255)
            sz = tmpframe.shape
            _, img_binary = cv2.threshold(tmpframe, 50, 255, cv2.THRESH_BINARY)
            _, contours, hierarchy = cv2.findContours(img_binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            contours = list(filter(lambda x: (cv2.contourArea(x) > 200)*(cv2.contourArea(x) < 650), contours))
            if len(contours) >= 1:
                cxy = []
                for i in range(len(contours)):
                    M = cv2.moments(contours[i])
                    cx = int(M['m10']/M['m00'])
                    cy = int(M['m01']/M['m00'])
                    cxy.append((cx,cy))
                if len(contours) == 1:
                    ps = (count,cxy[0][0],cxy[0][1])
                    ag = 0
                else:
                    cxy = np.array(cxy, dtype = float)
                    if np.sum(np.isnan(pos[-1])) == 0:
                        cxy_ = cxy - np.array(pos[-1][1:], dtype = float)
                    else:
                        cxy_ = cxy - np.array([sz[1],sz[0]], dtype = float)/2
                    ag = np.argmin(cxy_[:,0]**2 + cxy_[:,1]**2)
                    ps = (count,int(cxy[ag,0]),int(cxy[ag,1]))
                cv2.drawContours(frame_, contours[ag], -1, color=(0, 255, 0), thickness=3)
                cv2.circle(frame_, ps[1:], 10, (255,0,0), thickness=2, lineType=cv2.LINE_8, shift=0)  
                pos.append(ps)
            else:
            	pos.append((count,np.nan,np.nan))
            if param[6] > 0:
                if count == param[6]:
                    cv2.imshow("color",tmpframe)
                    cv2.waitKey(0)
                    return
            else:
                video.write(frame_)
            count += 1
            if count % 100 == 0:
            	print(count)
        else:
            print("video length",count)
            scipy.io.savemat(datapath + videoname2 + ".mat", {'pos':pos})
            video.release()
            cap.release()
            return

params = dict()
params[0] = (200,430,320,470,50,1.2,0)  #unable to detect
params[1] = (150,380,240,390,0,1.3,0)
params[2] = (100,330,80,230,100,1.2,0)
params[3] = (200,450,330,480,100,1.4,0)
params[4] = (130,400,330,520,140,1.25,0)  #short length
params[5] = (180,420,260,490,150,1.4,0)
params[6] = (80,320,360,550,150,1.4,0)
params[7] = (200,450,200,400,140,1.4,0)
params[8] = (70,370,200,400,130,1.4,0)
params[9] = (150,400,200,400,130,1.4,0)
params[10] = (150,400,300,500,130,1.3,0)   #short length?

subno = 7

patha = "/Users/TP/research/Tlab/behavior/human/20220502/"
pathb = "/Users/TP/research/Tlab/behavior/human/20220513/"

dr = os.listdir(patha)
dr2 = sorted([d for d in dr if d.endswith("avi")])
drb = os.listdir(pathb)
drb2 = sorted([d for d in drb if d.endswith("avi")])

dr2.extend(drb2)


if subno < 3:
    path = patha
else:
    path = pathb

targetv = dr2[subno]
print(targetv)

param = params[subno]

detect_object(path, targetv, targetv.replace("video","detect").split(".")[0], param)







