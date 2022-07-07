import cv2
import os

video_name = 'video.avi'

figpath = "/Users/yoshiki/research/Tlab/behavior/acceleration/20201207/movies/figures/"
vpath = "/Users/yoshiki/research/Tlab/behavior/acceleration/20201207/movies/"

dirs = os.listdir(figpath)
#rmv = [17,41,42,64,88,111,112,134,136,159,183,206,208,209,]
rmv = []
images = [[x for x in dirs if x.endswith("_"+str(i) + ".png")][0] for i in range(1,179) if i not in rmv]

image = cv2.imread(os.path.join(figpath, images[0]))
height, width, layers = image.shape

print(height, width)
fps = 24.75
fourcc = cv2.VideoWriter_fourcc(*'XVID')
video = cv2.VideoWriter(vpath + video_name, fourcc, fps, (567,height))


for image in images:
    img2 = cv2.imread(os.path.join(figpath, image))
    img2 = cv2.rectangle(img2, (54, 200), (63, 262), (0, 0, 0), -1)
    img2 = cv2.rotate(img2, cv2.ROTATE_90_CLOCKWISE)
    img2 = cv2.putText(img2, "2 m/s", (50,40), fontFace=cv2.FONT_ITALIC,
     fontScale = 0.8, color=(0,0,0),thickness=2,lineType=cv2.LINE_AA)
    img2 = cv2.putText(img2, "3", (137,33), fontFace=cv2.FONT_ITALIC,
     fontScale = 0.5, color=(0,0,0),thickness=2,lineType=cv2.LINE_AA)
    img2 = cv2.rotate(img2, cv2.ROTATE_90_COUNTERCLOCKWISE)
    img = cv2.resize(img2, (567, 319))
    video.write(img)

cv2.destroyAllWindows()
video.release()




