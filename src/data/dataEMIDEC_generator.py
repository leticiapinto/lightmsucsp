#leer carpeta de datasets
import os
from medpy.metric.binary import hd, dc
import nibabel as ni    
import matplotlib.pyplot as plot 
from skimage.io import imshow 

import cv2  # Not actually necessary if you just want to create an image.
import numpy as np
from pathlib import Path

def graph_generator(readfile, writefile):

    file_image = readfile + ".nii"
    
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

# importing the module
import os

# Providing the path of the directory
emidecpath = '/home/leticia/datasets/emidec-dataset-1.0.1'

# Retrieving the list of all the files
folders = os.listdir(emidecpath)

#folders = [a for a in os.listdir(emidecpath) if os.path.isdir(emidecpath)]

# loop to iterate every item in the list
'''
for file in folders:
    if Path(emidecpath+ "/" + file).is_dir():
        print(file)
'''
#for r, d, f in os.walk(emidecpath):

write_graph_dir = "/data/leticia/EMIDEC/graph_files/"

for f in folders:
        #print("f: ", f)
        if Path(emidecpath+ "/" + f).is_dir():
            print(f)
            #For images
            image_file          = emidecpath+ "/" + f + "/Images/" + f 
            groundt_file        = emidecpath+ "/" + f + "/Contours/" + f
            write_graph_file    = write_graph_dir + "gt_" +f
            print(image_file, groundt_file, write_graph_file)
            graph_generator(groundt_file, write_graph_file)


