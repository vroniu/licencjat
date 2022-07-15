# %%
import pandas as pd
import numpy as np
from keras.models import Sequential
from keras.metrics import Precision
from keras.layers import Dense, Conv2D, Flatten
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from converter import convert_x, convert_y, convert_y_predicted
from predicter_helper import convert_seq_to_model_input, prepare_primer3_input, parse_primer3_output
from random import randint, choice
import matplotlib.pyplot as plt
import os
# %%
dataset = pd.read_csv('dataset_from_orange.csv')
X = dataset.iloc[:, 0:20]
Y = dataset['score']
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, train_size=0.7)
X_train = convert_x(X_train)
X_test = convert_x(X_test)
Y_train = convert_y(Y_train)
Y_test = convert_y(Y_test)

model = Sequential([
    Conv2D(30, kernel_size=4, activation='relu', input_shape=(20, 4, 1)),
    Flatten(),
    Dense(30, activation='relu'),
    Dense(15, activation='relu'),
    Dense(5, activation='relu'),
    Dense(2, activation='softmax'),
])

model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
# %%
model.fit(X_train, Y_train, epochs=20)
# %%
Y_predicted = model.predict(X_test)
# %%
Y_predicted = convert_y_predicted(Y_predicted)
# %%
precision = Precision()
precision.update_state(y_pred=Y_predicted, y_true=Y_test)
print("using CNN - precision: " + str(precision.result().numpy()))

# with open('genome', 'r') as genome_input_file:
#     genome_string = ''.join(choice('ACTG') for _ in range(1100))
#     for random_attempt in range(10):
#         print('----------- PODEJŚCIE %d ------------' % random_attempt)
#         nt_count = 100
#         start_nt_index = randint(100, 1000)

#         with open('p3pick_input', 'w') as primer3_input_file:
#             primer3_input_file.write(prepare_primer3_input(genome_string, start_nt_index, nt_count))
#             primer3_input_file.close()
#         os.system('primer3_core --output=p3pick_output < p3pick_input')
#         with open('p3pick_output', 'r') as primer3_output_file:
#             twoja_stara = primer3_output_file.read()
#             nwm = parse_primer3_output(twoja_stara)
#             for entry in nwm:
#                 left_primer = entry[1]
#                 right_primer = entry[3]
#                 left_primer_pred = convert_y_predicted(model.predict(convert_seq_to_model_input(left_primer)))
#                 right_primer_pred = convert_y_predicted(model.predict(convert_seq_to_model_input(right_primer)))
#                 print('lewy primer ', end=" ")
#                 print(left_primer, end=" ")
#                 print('lewy primer predykcja ', end=" ")
#                 print(left_primer_pred)
#                 print('prawy primer ', end=" ")
#                 print(right_primer, end=" ")
#                 print('prawy primer predykcja ', end=" ")
#                 print(right_primer_pred)
#             primer3_output_file.close()
#     print('---------------------')

#     genome_input_file.close()
# %%
# Test
with open('genome', 'r') as genome_input_file:
    genome_string = genome_input_file.read().replace("\n", "")
    for random_attempt in range(9):
        print('----------- PODEJŚCIE %d ------------' % random_attempt)
        nt_count = 100
        start_nt_index = randint(100, len(genome_string) - 100)

        with open('p3pick_input', 'w') as primer3_input_file:
            primer3_input_file.write(prepare_primer3_input(genome_string, start_nt_index, nt_count))
            primer3_input_file.close()
        os.system('primer3_core --output=p3pick_output < p3pick_input')
        with open('p3pick_output', 'r') as primer3_output_file:
            primers_from_primer3 = parse_primer3_output(primer3_output_file.read())
            primer3_output_file.close()
        # search for left primers
        left_primer_end = start_nt_index - 1
        left_primer_found = []
        while left_primer_end > 20 and len(left_primer_found) < 5:
            potential_primer = genome_string[left_primer_end - 20 : left_primer_end]
            potential_primer_prediction = model.predict(convert_seq_to_model_input(potential_primer))
            if potential_primer_prediction[0][1] > 0.9:
                left_primer_found.append((potential_primer, potential_primer_prediction[0], left_primer_end - 20))
            left_primer_end = left_primer_end - 5
        # search for right primers
        right_primer_start = start_nt_index + nt_count
        right_primer_found = []
        while right_primer_start < len(genome_string) - 20 and len(right_primer_found) < 5:
            potential_primer = genome_string[right_primer_start : right_primer_start + 20]
            potential_primer_prediction = model.predict(convert_seq_to_model_input(potential_primer))
            if potential_primer_prediction[0][1] > 0.9:
                right_primer_found.append((potential_primer, potential_primer_prediction[0], right_primer_start))
            right_primer_start = right_primer_start + 5
        plotted_lines = 0.9
        plot_X = [x for x in range(start_nt_index, start_nt_index + nt_count)]
        plot_Y = [plotted_lines for x in plot_X]
        plt.plot(plot_X, plot_Y, color='grey', linewidth=10)
        plotted_lines = 1
        min_x, max_x = len(genome_string), 0
        for right_primer_start_pos in [x[2] for x in right_primer_found]:
            plot_X = [x for x in range(right_primer_start_pos, right_primer_start_pos + 20)]
            if plot_X[0] < min_x:
                min_x = plot_X[0]
            if plot_X[-1] > max_x:
                max_x = plot_X[-1]
            plot_Y = [plotted_lines for x in plot_X]
            plt.plot(plot_X, plot_Y, color='red', linewidth=10)
            plotted_lines = plotted_lines + 0.1
        plotted_lines = 1
        for left_primer_start_pos in [x[2] for x in left_primer_found]:
            plot_X = [x for x in range(left_primer_start_pos, left_primer_start_pos + 20)]
            if plot_X[0] < min_x:
                min_x = plot_X[0]
            if plot_X[-1] > max_x:
                max_x = plot_X[-1]
            plot_Y = [plotted_lines for x in plot_X]
            plt.plot(plot_X, plot_Y, color='red', linewidth=10)
            plotted_lines = plotted_lines + 0.1
        plotted_lines = 0.95
        for left_primer_start_pos in [int(x[0].split(",")[0]) for x in primers_from_primer3]:
            plot_X = [x for x in range(left_primer_start_pos, left_primer_start_pos + 20)]
            if plot_X[0] < min_x:
                min_x = plot_X[0]
            if plot_X[-1] > max_x:
                max_x = plot_X[-1]
            plot_Y = [plotted_lines for x in plot_X]
            plt.plot(plot_X, plot_Y, color='green', linewidth=10)
            plotted_lines = plotted_lines + 0.1
        plotted_lines = 0.95
        for right_primer_start_pos in [int(x[2].split(",")[0]) for x in primers_from_primer3]:
            plot_X = [x for x in range(right_primer_start_pos, right_primer_start_pos + 20)]
            if plot_X[0] < min_x:
                min_x = plot_X[0]
            if plot_X[-1] > max_x:
                max_x = plot_X[-1]
            plot_Y = [plotted_lines for x in plot_X]
            plt.plot(plot_X, plot_Y, color='green', linewidth=10)
            plotted_lines = plotted_lines + 0.1
        plt.xlim([min_x - 20, max_x + 20])
        ax = plt.gca()
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        plt.show()
        print('---------------------')

    genome_input_file.close()
