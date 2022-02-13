from operator import le
from statistics import mean, stdev
import numpy
from sklearn import preprocessing
from sklearn.model_selection import StratifiedKFold
from sklearn import linear_model
from sklearn import datasets
from data_csv_utils import read_data, convert
import numpy as np
from KNN.main import train as knn_train
from encoding import one_hot_encode_sequences, mapping_encode_value

def shuffle_and_train():
    kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=12)
    location, status, fragment = read_data()
    # encode status and fragment
    encoded_status = []
    encoded_fragment = []
    encode_mapper, decode_mapper = one_hot_encode_sequences()
    for i, str in enumerate(fragment):
        temp = []
        for j, c in enumerate(str):
            temp.append(mapping_encode_value(encoding_mapper=encode_mapper, decoding_mapper=decode_mapper, key=c))
        encoded_fragment.append(np.array(temp).flatten())
    for i, str in enumerate(status):
        if(str == 'nonSNO'):
            encoded_status.append(0)
        else:
            encoded_status.append(1)
    data_train = []
    label_train = []
    data_test = []
    label_test = []
    for train_index, test_index in kf.split(encoded_fragment, encoded_status):
        data_train.append(np.array(encoded_fragment)[train_index])
        label_train.append(np.array(encoded_status)[train_index])
        data_test.append(np.array(encoded_fragment)[test_index])
        label_test.append(np.array(encoded_status)[test_index])
        # print(data_train.shape)
        # print(label_train)
        # print(data_test.shape)
        # print(np.array(encoded_status).shape)
    knn_train(data_train, label_train, data_test, label_test, 5)
    #     train_one_fold(data_train, label_train)
    #     test_one_fold(data_test, label_test)

shuffle_and_train()
