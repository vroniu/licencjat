# %%
# Using word embedding
import pandas as pd
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Flatten
from keras.layers.embeddings import Embedding
from keras.metrics import Precision

# %%
dataset = pd.read_csv('dataset_from_orange.csv')
for i in range(20):
    dataset['data_%d' % (i + 1)] = dataset['data_%d' % (i + 1)].astype(int)
# %%
X_train, X_test, Y_train, Y_test = train_test_split(
    dataset.iloc[:, 0:20],
    dataset.iloc[:, 20],
    train_size=0.7
)
# %%
model = Sequential([
    Embedding(5, 80, input_length=20),
    Flatten(),
    Dense(60),
    Dense(40),
    Dense(20),
    Dense(10),
    Dense(1, activation='softmax')
])
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy', 'precision'])
# %%
# history = model.fit(X_train, Y_train, epochs=100)
model.load_weights('word_embedding_weights')
# %%
Y_pred = model.predict(X_test)
precision = Precision()
precision.update_state(y_pred=Y_pred, y_true=Y_test)
print("Using word embedding - precision: " + str(precision.result().numpy()))
# %%
model.save_weights('word_embedding_weights')
# %%
