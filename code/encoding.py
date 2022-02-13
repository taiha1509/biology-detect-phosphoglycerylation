from base64 import decode
from sklearn.preprocessing import OneHotEncoder


def one_hot_encode_sequences(
        X = [["V"],["X"],["Y"],
        ["D"],["T"],["S"],["L"],["I"],
        ["W"],["H"],["Q"],["N"],["G"],
        ["A"],["K"],["M"],["R"],["P"],
        ["C"],["F"],["E"], ["-"]]
    ):
    x_fit = OneHotEncoder().fit(X)
    encode = x_fit.transform(X)
    return encode.toarray(), X

def mapping_encode_value(encoding_mapper, decoding_mapper, key):
    index = 0
    for i, item in enumerate(decoding_mapper):
        if(item[0] == key):
            index = i
            break
    return encoding_mapper[index]
