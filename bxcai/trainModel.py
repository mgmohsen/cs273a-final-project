import clean as c
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn.metrics import roc_auc_score

import numpy as np

k = 6
regionList = ['positive', 'heart', 'neural tube', 'limb', 'brain', 'hindbrain', 'midbrain', 'forebrain']

rawData = c.parseRefData(k)
negativeRatio = 1

for region in ['limb']:
	X, y = c.processSVMData(rawData, region)
	pi = [i for i, v in enumerate(y) if v == 1]
	ni = [i for i, v in enumerate(y) if v == 0]
	ns = np.random.choice(ni, size=len(pi)*negativeRatio)
	indices = np.concatenate((pi, ns))
	X = X[indices]
	y = y[indices]
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5)
	clf = svm.SVC(kernel='linear', probability=True)
	clf.fit(X_train, y_train)
	print(region)
	print("train roc: {}".format(roc_auc_score(y_train, clf.decision_function(X_train))))
	print("test roc: {}".format(roc_auc_score(y_test, clf.decision_function(X_test))))
