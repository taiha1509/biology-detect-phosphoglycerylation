from cProfile import label
from os import stat
import scipy.io as sio
import csv

def read_data():
    f_data = open('../data/data.csv', 'r')
    csv_reader = csv.reader(f_data, delimiter = ',')
    location = []
    status = []
    fragment = []
    index = 0
    for row in csv_reader:
        if (len(row) > 0):
            if(index == 0):
                index += 1
                continue
            else:
                index +=1
                location.append(row[0])
                status.append(row[1])
                fragment.append(row[2])

    return location, status, fragment

def convert():
    f_data = open('./data/data.csv', 'w')
    writer = csv.writer(f_data)
    header = ['Location', 'Status', 'Fragment']
    writer.writerow(header)
    location = []
    status = []
    fragment = []

    file_dataset = "../RAM-PGK-master/Phosphoglycerylationstruct.mat"
    dataset = sio.loadmat(file_dataset)

    all_sequences = dataset["DB_Phosphoglycerylation"][0]
    gaps = '--------------------'
    temp_sequence = ''
    for i, s in enumerate(all_sequences):
        temp_sequence = gaps + s[0][0][0][0] + gaps
        for j, amino_axit in enumerate(temp_sequence):
            if(amino_axit == 'K'):
                location.append(j - 19)
                label = s[3][0][0][0][j-20]
                if(label == '1'):
                    status.append('SNO')
                else:
                    status.append('nonSNO')
                fragment.append(temp_sequence[j-20:j+21])

    for i, c in enumerate(location):
        writer.writerow([location[i], status[i], fragment[i]])
    f_data.close()
