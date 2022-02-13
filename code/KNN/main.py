import enum
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold
import numpy as np
from sklearn import linear_model, datasets, preprocessing
import sys
sys.path.append('../')
from data_csv_utils import read_data, convert
from encoding import one_hot_encode_sequences, mapping_encode_value
import math

def divide_dataset():
    location, status, fragment = read_data()
    
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
    X_train, X_test, y_train, y_test = train_test_split(encoded_fragment, encoded_status, test_size=0.33, random_state=42, shuffle=True)
    return X_train, X_test, y_train, y_test

def calc_average(lst):
    if(len(lst == 0)):
        return 0
    return sum(lst) / len(lst)

def train(X_train, label_train, X_validation, label_validation, K_Fold):
    neigh_model_list = []
    for i in range(0, K_Fold - 1):
        neigh = KNeighborsClassifier(n_neighbors=1)
        neigh.fit(X_train[i], label_train[i])
        Y = neigh.predict(X_validation[i])
        # danh gia mo hinh
        Accuracy_dict = {'Sensitivity': [], 'Specificity': [], 'Precision': [], 'Accuracy': [], 'MCC': []}
        TP = 0; FP = 0; FN = 0; TN = 0
        for j, y in enumerate(Y):
            if(y == 1): 
                if(label_validation[i][j] == 1):
                    TP += 1
                else:
                    FP += 1
            else:
                if(label_validation[i][j] == 1):
                    FN += 1
                else:
                    TN += 1
        # Sensitivity
        if((TP + FN) != 0):
            Accuracy_dict['Sensitivity'].append(TP/(TP + FN))
        else:
            Accuracy_dict['Sensitivity'].append(0)
        #Specificity
        if((TN + FP) != 0):
            Accuracy_dict['Specificity'].append(TN/(TN + FP))
        else:
            Accuracy_dict['Sensitivity'].append(0)
        #Precision
        if((TP + FP) != 0):
            Accuracy_dict['Precision'].append(TP/(TP + FP))
        else:
            Accuracy_dict['Sensitivity'].append(0)
        #Accuracy
        Accuracy_dict['Accuracy'].append((TN + TP)/(TN + FP + FN + TP))
        #MCC
        if((TP+FP)*(TP+FN)*(TN+FP)*(TN+FP) != 0):
            Accuracy_dict['MCC'].append(((TN * TP) - (FN * FP))/math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FP)))
        else:
            Accuracy_dict['MCC'].append(0)
        neigh_model_list.append(neigh)

    print("Sensitivity: ", np.average(Accuracy_dict['Sensitivity']))
    print("Specificity: ", np.average(Accuracy_dict['Specificity']))
    print("Precision: ", np.average(Accuracy_dict['Precision']))
    print("Accuracy: ", np.average(Accuracy_dict['Accuracy']))
    print("MCC: ", np.average(Accuracy_dict['MCC']))
        
    return neigh_model_list

def test(model_list, X_test, y_test, K_Fold):
    print('###########################################', np.array(X_test).shape)
    output_label = []
    label_after_voting = []
    for i, model in enumerate(model_list):
        Y = model.predict(X_test)
        output_label.append(Y)
    for i, label_of_first_model in enumerate(output_label[0]):
        nonSNO_voting_num = 0
        for j in range (0, K_Fold - 1):
            if(output_label[j][i] == 0):
                nonSNO_voting_num += 1
        if(nonSNO_voting_num >= 3):
            label_after_voting.append(1)
        else:
            label_after_voting.append(0)


    # danh gia mo hinh
    Accuracy_dict = {'Sensitivity': [], 'Specificity': [], 'Precision': [], 'Accuracy': [], 'MCC': []}
    TP = 0; FP = 0; FN = 0; TN = 0
    for j, y in enumerate(label_after_voting):
        if(y == 1): 
            if(y_test[j] == 1):
                TP += 1
            else:
                FP += 1
        else:
            if(y_test[j] == 1):
                FN += 1
            else:
                TN += 1
    # Sensitivity
    if((TP + FN) != 0):
        Accuracy_dict['Sensitivity'].append(TP/(TP + FN))
    else:
        Accuracy_dict['Sensitivity'].append(0)
    #Specificity
    if((TN + FP) != 0):
        Accuracy_dict['Specificity'].append(TN/(TN + FP))
    else:
        Accuracy_dict['Sensitivity'].append(0)
    #Precision
    if((TP + FP) != 0):
        Accuracy_dict['Precision'].append(TP/(TP + FP))
    else:
        Accuracy_dict['Sensitivity'].append(0)
    #Accuracy
    Accuracy_dict['Accuracy'].append((TN + TP)/(TN + FP + FN + TP))
    #MCC
    if((TP+FP)*(TP+FN)*(TN+FP)*(TN+FP) != 0):
        Accuracy_dict['MCC'].append(((TN * TP) - (FN * FP))/math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FP)))
    else:
        Accuracy_dict['MCC'].append(0)

    print("Sensitivity: ", np.average(Accuracy_dict['Sensitivity']))
    print("Specificity: ", np.average(Accuracy_dict['Specificity']))
    print("Precision: ", np.average(Accuracy_dict['Precision']))
    print("Accuracy: ", np.average(Accuracy_dict['Accuracy']))
    print("MCC: ", np.average(Accuracy_dict['MCC']))


def shuffle_and_train(K_Fold = 5):
    X_train, X_test, y_train, y_test = divide_dataset()
    kf = StratifiedKFold(n_splits=K_Fold, shuffle=True, random_state=12)
    data_train = []
    label_train = []
    data_validation = []
    label_validation = []
    for train_index, test_index in kf.split(X_train, y_train):
        data_train.append(np.array(X_train)[train_index])
        label_train.append(np.array(y_train)[train_index])
        data_validation.append(np.array(X_train)[test_index])
        label_validation.append(np.array(y_train)[test_index])
        # print(data_train.shape)
        # print(label_train)
        # print(data_test.shape)
        # print(np.array(encoded_status).shape)
    neigh_model_list = train(data_train, label_train, data_validation, label_validation, K_Fold)
    test(neigh_model_list, X_test, y_test, K_Fold)

shuffle_and_train()
