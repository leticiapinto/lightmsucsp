import os 
import glob
from pathlib import Path


def get_train_val_paths(input_folder, split_train_val = True):

    file_list = {'val': [], 'train': []}
    for folder in os.listdir(input_folder):
        folder_path = os.path.join(input_folder, folder)
        if os.path.isdir(folder_path):
            if split_train_val:
                train_test = 'val' if (int(folder[-3:]) % 5 == 0) else 'train'
            else:
                train_test = 'train'
            search = 'patient???_frame??.nii.gz'

            for file in glob.glob(os.path.join(folder_path, search)):
                file_list[train_test].append(Path(file))
    return file_list

def get_da_files(dadirectory):
    lstda = os.listdir(dadirectory)
    return lstda


def get_da_bbox(filename):

    f = open(filename, "r")

    lst_da_bbox = []
    for line in f:
        #print(line.split("\t"))
        line_bbox = line.split("\t")
        #### read by rows and check if it doesn't exist then save it
        min_x = line_bbox[5]
        min_y = line_bbox[6]
        max_x = line_bbox[7]
        max_y = line_bbox[8].rstrip()
        exists = False

        for i in range(0, len(lst_da_bbox)):
            if(lst_da_bbox[i][0] == min_x and lst_da_bbox[i][1] == min_y and lst_da_bbox[i][2] == max_x and lst_da_bbox[i][3] == max_y):
                ## This one already exists
                exists = True 
                continue
        
        if(not exists):
            lst_bbox = []
            #lst_bbox.append(line_bbox[0])
            #lst_bbox.append(line_bbox[1])
            lst_bbox.append(min_x)
            lst_bbox.append(min_y)
            lst_bbox.append(max_x)
            lst_bbox.append(max_y)
            lst_da_bbox.append(lst_bbox)

    return lst_da_bbox
    