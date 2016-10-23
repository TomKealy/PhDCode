import pandas as pd
import numpy as np
from sklearn import cross_validation
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
from sklearn.linear_model import SGDClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
from nolearn.dbn import DBN
import timeit

train = pd.read_csv("mnist_train.csv")
features = train.columns[1:]
X = train[features]
y = train['label']
X_train, X_test, y_train, y_test = cross_validation.train_test_split(X/255., y, test_size = 0.1, random_state = 0)

clf_rf = RandomForestClassifier()
clf_rf.fit(X_train, y_train)

y_pred_rf = clf_rf.predict(X_test)
acc_rf = accuracy_score(y_test, y_pred_rf)

"nn"

clf_nn = DBN([X_train.shape[1], 300, 10], learn_rates=0.3, learn_rate_decays=0.9, epochs=50)

clf_nn.fit(X_train, y_train)
acc_nn = clf_nn.score(X_test, y_test)


