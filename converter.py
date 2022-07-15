import numpy as np
from sklearn.preprocessing import OneHotEncoder

def convert_x(X):
    X_2d = []
    for primer_data in X.values:
        primer_data_2d = []
        for i in range(20):
            primer_data_at_current_position = [0, 0, 0, 0]
            nucleotide = primer_data[i]
            primer_data_at_current_position[int(nucleotide) - 1] = 255
            primer_data_2d.append(primer_data_at_current_position)
        X_2d.append(primer_data_2d)
    X_2d = np.array(X_2d)
    X_2d = X_2d.reshape(len(X), 20, 4, 1)
    return X_2d

def convert_y(Y):
    encoder = OneHotEncoder(sparse=False)
    Y = np.array(Y)
    Y = Y.reshape(-1, 1)
    Y_2d = encoder.fit_transform(Y)
    return Y_2d

def convert_y_predicted(Y_predicted):
    Y_predicted_roundedup = []
    for prediction in Y_predicted:
        if prediction[0] > prediction[1]:
            Y_predicted_roundedup.append([1, 0])
        else:
            Y_predicted_roundedup.append([0, 1])
    Y_predicted = np.array(Y_predicted_roundedup)
    return Y_predicted_roundedup