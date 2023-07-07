from PIL import Image
import numpy as np
import IPython.display as fig
import matplotlib.pyplot as plt
from torchvision.utils import make_grid, save_image
import torch
import IPython.display as fig
from io import BytesIO
import time, os
from torchvision import transforms, utils, models
import cv2
from torchvision import transforms, utils


import os, json
import cv2
import torch.utils.data
from PIL import Image
import numpy as np
import re
import glob
from pathlib import Path
from torchvision import transforms, utils
import nibabel as nib
from collections import namedtuple
import reading 


def generatepatientlist():
    #read patientsfile.txt
    num_files    = 1900
    filename = 'patientsfile.txt'
    f = open(filename, "r")
    with open(dataaugmentation + patientfile,) as pfile:
        lines = pfile.readlines()
    lst = [line.split("\n")[0] for line in lines]
    # 5, 50, 100, 500, 1000
    lstfiles = []
    lstfiles.append(lst[5])
    lstfiles.append(lst[50])
    lstfiles.append(lst[100])
    lstfiles.append(lst[500])
    lstfiles.append(lst[1000])

    return lstfiles


def generateimageandmask(list_images, folder_path):

    #images = os.listdir(folder_path  + 'image/' )
    ## I need to open NIFTI files
    
    #img_files = [folder_path + 'image/' +  x for x in images]
    trans = transforms.ToTensor()
    #/data/leticia/training/patient090/patient090_frame04.nii.gz
    #for index in range(0, 5):
    #for index in range(0, len(img_files)):
    index = 3
    count = 0 
    csize       = 128
    transform_ = transforms.Compose([
                                        transforms.Resize([int(csize),int(csize)]),
                                        transforms.ToTensor()
                                        ])
    for image_path in list_images:
        count += 1
        ##loading images and masks
        image_info = nib.load(image_path)
        image_ = np.array(image_info.dataobj)[:,:,index]
        cv2.imwrite(folder_path + "image_" + str(count) + ".jpg", image_) 

        imgname  = str(image_path).split(".")[0]
        mask_path = imgname + "_gt.nii.gz"
        mask_info = nib.load(mask_path)
        mask_ = np.array(mask_info.dataobj)[:,:,index]
        mask_ = np.where(mask_ == 3, 255, 0.0)
        #mask_ = transform_(Image.fromarray(mask_))
        cv2.imwrite(folder_path + "mask_" + str(count) + ".jpg", mask_)


def visualizeallimagesHGBUNet(list_images, model, folder_path, model_name):

    model.cuda()
    #images = os.listdir(folder_path  + 'image/' )
    ## I need to open NIFTI files
    
    #img_files = [folder_path + 'image/' +  x for x in images]
    trans = transforms.ToTensor()
    #/data/leticia/training/patient090/patient090_frame04.nii.gz
    #for index in range(0, 5):
    #for index in range(0, len(img_files)):
    index = 3
    count = 0 
    csize       = 128
    transform_ = transforms.Compose([
                                        transforms.Resize([int(csize),int(csize)]),
                                        transforms.ToTensor()
                                        ])
    for image_path in list_images:
        count += 1
        ##loading images and masks
        image_info = nib.load(image_path)
        image_ = np.array(image_info.dataobj)[:,:,index]
        #cv2.imwrite(folder_path + "image_" + str(count) + ".jpg", image_) 

        #imgname  = str(image_path).split(".")[0]
        #mask_path = imgname + "_gt.nii.gz"
        #mask_info = nib.load(mask_path)
        #mask_ = np.array(mask_info.dataobj)[:,:,index]
        #mask_ = np.where(mask_ == 3, 255, 0.0)
        #mask_ = transform_(Image.fromarray(mask_))
        #cv2.imwrite(folder_path + "mask_" + str(count) + ".jpg", mask_) 

        ##applying preprocess
        max_image = np.amax(image_)
        image_ = (image_/max_image)
        #print(np.amax(image))
        image_ = transform_(Image.fromarray(image_))
        temp = torch.zeros(3, 128, 128, dtype=torch.float)
        temp[0] = image_
        temp[1] = image_
        temp[2] = image_
        image_ = temp
        #print(image.shape)
        #print(np.unique(mask))
        

        output  = model(image_.unsqueeze(0).cuda()).cpu()
        output  = torch.sigmoid(output).detach()[0].numpy()
        output[output < 0.5] = 0
        output[output >= 0.5] = 255
        #print('len(output > 1.0) ', len(output[output > 1.0 ]))
        #print('len(output < 0) ', len(output[output < 0 ]))
        #print('output.shape ', output.shape)
        output = np.transpose(output,  (1, 2, 0))
        #print('output.shape ', output.shape)
        image_ = image_.cpu().detach().numpy()
        #print(type(image_), image_.shape)
        #print(type(image_), image_.shape)
        
        cv2.imwrite(folder_path + model_name + "_" + str(count) + ".jpg", output) 

def visualizeallimagesHGB_FCNvgg(list_images, model, folder_path, model_name):

    model.cuda()
    index = 3
    count = 0 
    csize       = 128
    transform_ = transforms.Compose([
                                        transforms.Resize([int(csize),int(csize)]),
                                        transforms.ToTensor()
                                        ])

    for image_path in list_images:
        count += 1
        image_info = nib.load(image_path)
        image_ = np.array(image_info.dataobj)[:,:,index]

        max_image = np.amax(image_)
        image_ = (image_/max_image)
        #print(np.amax(image))
        image_ = transform_(Image.fromarray(image_))
        temp = torch.zeros(3, 128, 128, dtype=torch.float)
        temp[0] = image_
        temp[1] = image_
        temp[2] = image_
        image_ = temp

        #print(img_.unsqueeze(0).shape)
        output  = model(image_.unsqueeze(0).cuda()).cpu()
        output  = torch.sigmoid(output).detach()[0].numpy()
        #print('output.shape ', output.shape)
        output[output < 0.5] = 0
        output[output >= 0.5] = 255
        #output = np.transpose(output,  (1, 2, 0))
        #print('output.shape ', output.shape)
        #cv2.imwrite('results/'+ dataname + '/' + images[index], output)
        cv2.imwrite(folder_path + model_name + "_" + str(count) + ".jpg", output[0]) 

def visualizeallimagesHGB_FCNResNet(list_images, model, folder_path, model_name):

    model.cuda()
    index = 3
    count = 0 
    csize       = 128
    transform_ = transforms.Compose([
                                        transforms.Resize([int(csize),int(csize)]),
                                        transforms.ToTensor()
                                        ])

    for image_path in list_images:
        count += 1
        image_info = nib.load(image_path)
        image_ = np.array(image_info.dataobj)[:,:,index]

        max_image = np.amax(image_)
        image_ = (image_/max_image)
        #print(np.amax(image))
        image_ = transform_(Image.fromarray(image_))
        temp = torch.zeros(3, 128, 128, dtype=torch.float)
        temp[0] = image_
        temp[1] = image_
        temp[2] = image_
        image_ = temp
        
        output  = model(image_.unsqueeze(0).cuda())['out'].cpu()
        output  = torch.sigmoid(output).detach()[0].numpy()
        output  = np.asarray(output)[0]
        output[output < 0.5] = 0
        output[output >= 0.5] = 255
        print(type(output) , output.shape)

        #output = np.transpose(output,  (1, 2, 0))
        #print('output.shape ', output.shape)
        #image_ = image_.cpu().detach().numpy()

        cv2.imwrite(folder_path + model_name + "_" + str(count) + ".jpg", output)


def plots(dataname, val_losses, train_losses, val_iou, train_iou, val_dice, train_dice, val_prec, train_prec, val_jc, train_jc):

    plt.figure(figsize = (24, 5))
    plt.subplot(1, 5, 1)
    plt.plot(val_losses, 'bo-', label = 'val-loss')
    plt.plot(train_losses, 'ro-', label = 'train-loss')
    plt.grid('on')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['validation', 'training'], loc='upper left')

    plt.subplot(1, 5, 2)
    plt.plot(val_iou, 'bo-', label = 'val-iou')
    plt.plot(train_iou, 'ro-', label = 'train-iou')
    plt.ylabel('Intersection Over Union')
    plt.grid('on')
    plt.xlabel('epoch')
    plt.legend(['validation', 'training'], loc='upper left')

    plt.subplot(1, 5, 3)
    plt.plot(val_dice, 'bo-', label = 'val-dice')
    plt.plot(train_dice, 'ro-', label = 'train-dice')
    plt.ylabel('Dice')
    plt.grid('on')
    plt.xlabel('epoch')
    plt.legend(['validation', 'training'], loc='upper right')

    plt.subplot(1, 5, 4)
    plt.plot(val_prec, 'bo-', label = 'val-prec')
    plt.plot(train_prec, 'ro-', label = 'train-prec')
    plt.ylabel('Precision')
    plt.grid('on')
    plt.xlabel('epoch')
    plt.legend(['validation', 'training'], loc='upper right')

    plt.subplot(1, 5, 5)
    plt.plot(val_jc, 'bo-', label = 'val-jc')
    plt.plot(train_jc, 'ro-', label = 'train-jc')
    plt.ylabel('Jaccard')
    plt.grid('on')
    plt.xlabel('epoch')
    plt.legend(['validation', 'training'], loc='upper right')
    plt.savefig( dataname + '_plots.png')
    plt.show()