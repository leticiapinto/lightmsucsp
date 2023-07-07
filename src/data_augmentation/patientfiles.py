import os
num_lines = 0   #### this is the number of lines or bboxes ...
num_files = 5
path = '/data/leticia/training_da/'

file = open('patientsfile.txt', 'a')

for r, d, f in os.walk(path):
    for file_ in f:
        ##read files and count lines
        print(file_)
        #a,b = file.split(".")
        file.write(file_ + '\n')

file.close()