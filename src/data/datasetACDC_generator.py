#leer carpeta de datasets
import os
from medpy.metric.binary import hd, dc
import nibabel as ni    
import matplotlib.pyplot as plot 
from skimage.io import imshow 

import cv2  # Not actually necessary if you just want to create an image.
import numpy as np

def graph_generator(readfile, writefile):

    file_image = readfile + ".nii.gz"
    
    #print(file_image,final_file)
    img = ni.load(file_image)
    img = img.get_data()
    s = img.shape
    rs = s[0]
    cs = s[1]
    max_range = s[2]
    #print(s)
    for k in range(max_range):
        final_file = writefile + "_" + str(k) +  ".graph"
        file = open(final_file,"w") 
        
        #print(max)     
        stline = '#rs ' + str(s[1]) + ' cs ' + str(s[0])+ '\n'
        file.write(stline)
        #print(stline)
        vertices = s[0]*s[1]
        edges = (s[0]-1)*( (s[1]-1)*2 + 1) + s[1] -1
        ndline = str(vertices) + ' ' + str(edges) + '\n'
        #print(ndline)
        file.write(ndline)
        message = 'val' + ' ' + 'sommets' + '\n'
        file.write(message)
        
        for i in range (vertices):
            ms = str(i) + ' ' + str('1') + '\n'
            file.write(ms)
            #print(ms)

        arcs = 'arcs values'
        file.write(arcs)
        #print(arcs)

        for i in range (rs-1):
            for j in range (cs-1):
                #print(img[i][j][k])
                tmp_h = str(i*cs + j+1)+ ' ' + str(i*cs + j)  + ' ' + str( abs(img[i][j+1][k] - img[i][j][k]))
                tmp_v = str((i+1)*cs + j) + ' ' + str(i*cs + j) + ' ' + str( abs(img[i+1][j][k] - img[i][j][k]))
                tmp = '\n' + tmp_h + '\n' + tmp_v
                #print(tmp)
                file.write(tmp)
            j = cs -1    
            tmp_v = '\n' + str((i+1)*cs + j) + ' ' + str(i*cs + j) + ' ' + str( abs(img[i+1][j][k] - img[i][j][k]))
            #print(tmp_v)
            file.write(tmp_v)

        i = rs - 1
        for j in range (cs-1):
            tmp_h = '\n' + str(i*cs + j+1) + ' ' + str(i*cs + j) + ' ' + str( abs(img[i][j+1][k] - img[i][j][k]))
            #print(tmp_h)
            file.write(tmp_h)

        file.close()

def image_generator(readfile, writefile):

    file_image = readfile + ".nii.gz"

    img = ni.load(file_image)
    img = img.get_data()
    s = img.shape
    rs = s[0]
    cs = s[1]
    max_range = s[2]
    
    for k in range(max_range):
        final_file = writefile + "_" + str(k) +  ".png"
        
        writeimage = np.zeros((rs,cs,1), np.uint8)
        writeimage = img[:,:,k]
        cv2.imwrite(final_file, writeimage)
    #print(file_image,final_file)


def gt_image_generator(readfile, writefile):

    file_image = readfile + ".nii.gz"

    img = ni.load(file_image)
    img = img.get_data()
    s = img.shape
    rs = s[0]
    cs = s[1]
    max_range = s[2]
    print(rs)
    print(cs)
    print(max_range)
    
    for k in range(max_range):
        writeimage = np.zeros((rs,cs,1), np.uint8)
        for i in range (rs-1):
            for j in range (cs-1):
                    writeimage[i][j] = img[i][j][k]*85
        final_file = writefile + "_" + str(k) +  ".png"
        cv2.imwrite(final_file, writeimage)

def gt_general_image_generator(readfile, writefile):

    file_image = readfile + ".nii.gz"

    img = ni.load(file_image)
    img = img.get_data()
    s = img.shape
    rs = s[0]
    cs = s[1]
    max_range = s[2]
    print(rs)
    print(cs)
    print(max_range)
    
    for k in range(max_range):
        writeimage = np.zeros((rs,cs,1), np.uint8)
        for i in range (rs-1):
            for j in range (cs-1):
                    writeimage[i][j] = img[i][j][k]
        final_file = writefile + "_" + str(k) +  ".png"
        cv2.imwrite(final_file, writeimage)        

#leer carpeta de datasets
# Getting the current work directory (cwd)
# r=root, d=directories, f = files
#os.makedirs("res/graph_files")
#os.makedirs("res/image_files")
training_dir = "training/patient"
write_graph_dir = "res/graph_files"
write_image_dir = "res/image_files"
write_gt_dir = "res/gt_files"

for patient in range(1,101):
    if patient < 10:
        patient = str('00') + str(patient)
    elif patient < 100:
        patient = str('0') + str(patient)
    #print(training_dir+str(k))    
    for r, d, f in os.walk(training_dir+str(patient)):
        for file in f:
            if ".nii" and "frame" in file:
                a,b,c = file.split(".")
                read_file = r + "/" + a
                write_graph_file = write_graph_dir + "/" + a
                write_image_file = write_image_dir + "/" + a
                write_gt_file = write_gt_dir + "/" + a

                if a.endswith("_gt"):
                    #usando el proceso normalizado para los tipo ground truth
                    print (a)
                    #gt_image_generator(read_file, write_gt_file)
                    gt_general_image_generator(read_file, write_gt_file) #mantiene etiquetas 0,1,2 y 3
                else:
                    #usando el proceso normal de generacion de graphs y de files
                    print("don't have gt: " + a)
                    #graph_generator(read_file, write_graph_file)
                    #image_generator(read_file, write_image_file) #para generar imagenes en grayscale
                a = a + '\n'
