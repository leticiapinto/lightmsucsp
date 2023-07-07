## read each file and count files to know how many lines does it have and how many different per file
dataaugmentation = '/data/leticia/training_da/'
filename = 'patientsfile.txt'

f = open(filename, "r")
num_patients = 0
num_files    = 1900
lst_da_bbox = []
for line in f:
    #print(line.split("\t"))
    patientfile = line.split("\n")
    patientfile = patientfile[0]
    #print(patientfile)
    with open(dataaugmentation + patientfile,) as pfile:
        lines = pfile.readlines()
        num_patients += len(lines)

print('total number of patients: ', num_patients)
print(num_patients/num_files)
    
