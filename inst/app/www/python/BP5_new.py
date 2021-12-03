# -*- coding: utf-8 -*-
from numpy.random import RandomState
from sklearn.decomposition import pca
from sklearn import model_selection
from sklearn.preprocessing import scale
import argparse
import numpy as np
import pandas as pd
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

# PCA
def get_pca(K, data):
    model = pca.PCA(n_components=int(K)).fit(data)
    data = model.transform(data)
    return data

# Turn labels into dataframe
def one_hot_dataframe(labels):
    #labels = labels.tolist()
    dict_idx = sorted(set(labels))
    one_hot_df = pd.DataFrame(columns=dict_idx)
    for i in range(len(labels)):
        one_hot_df.loc[i, labels[i]] = 1
    one_hot_df = one_hot_df.replace(to_replace=np.nan, value=0)
    return one_hot_df, dict_idx

def get_BP5_res(X_total_path, Y_total, X_verify, subclusters, num_classes, input_size, hidden_units_size1, hidden_units_size2, hidden_units_size3, regularizer1, regularizer2, learn_rate):
    num_classes = int(num_classes)
    input_size = int(input_size)
    hidden_units_size1 = int(hidden_units_size1)
    hidden_units_size2 = int(hidden_units_size2)
    hidden_units_size3 = int(hidden_units_size3)
    batch_size = 16
    training_iterations = 10000
    regularizer1 = regularizer1
    regularizer2 = regularizer2
    learn_rate = learn_rate


    X = tf.placeholder(tf.float32, shape=(None, input_size), name='x-input')
    Y = tf.placeholder(tf.float32, shape=(None, num_classes), name='y-input')

    W1 = tf.Variable(tf.random_normal([input_size, hidden_units_size1], stddev=0.1))
    B1 = tf.Variable(tf.constant(0.01), [hidden_units_size1])

    W2 = tf.Variable(tf.random_normal([hidden_units_size1, hidden_units_size2], stddev=0.1))
    B2 = tf.Variable(tf.constant(0.01), [hidden_units_size2])

    W3 = tf.Variable(tf.random_normal([hidden_units_size2, hidden_units_size3], stddev=0.1))
    B3 = tf.Variable(tf.constant(0.01), [hidden_units_size3])

    W4 = tf.Variable(tf.random_normal([hidden_units_size3, num_classes], stddev=0.1))
    B4 = tf.Variable(tf.constant(0.01), [num_classes])

    hidden_opt1 = tf.matmul(X, W1) + B1
    hidden_opt1 = tf.nn.relu(hidden_opt1)
    hidden_opt2 = tf.matmul(hidden_opt1, W2) + B2
    hidden_opt2 = tf.nn.relu(hidden_opt2)
    hidden_opt3 = tf.matmul(hidden_opt2, W3) + B3
    hidden_opt3 = tf.nn.relu(hidden_opt3)
    final_opt = tf.matmul(hidden_opt3, W4) + B4

    y = tf.nn.softmax(final_opt)

    loss = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=Y, logits=final_opt))+regularizer1*tf.nn.l2_loss(W1)+regularizer1*tf.nn.l2_loss(W2)+regularizer2*tf.nn.l2_loss(W3)+regularizer2*tf.nn.l2_loss(W4)

    opt = tf.train.AdamOptimizer(learn_rate).minimize(loss)

    correct_prediction = tf.equal(tf.argmax(Y, 1), tf.argmax(final_opt, 1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, 'float'))

    X_total =  pd.read_hdf(X_total_path, key='dge')
    X_total = X_total.T
    print(X_total.columns)
    #X_total = X_total.drop(['Unnamed: 0'], axis=1)

    Y_total, dict_idx = one_hot_dataframe(Y_total)
    print(Y_total.columns)
    print("get data")

    X_verify = X_verify
    print(X_verify.columns)
    #X_verify = X_verify.drop(['Unnamed: 0'], axis=1)
    a = X_verify.shape
    print()

    X_tmp = pd.concat([X_verify, X_total])
    X_tmp = X_tmp.replace(np.nan, 0)
    print(X_tmp.shape)
    #X_tmp = X_tmp.loc[:,~(X_tmp==0).all(axis=0)]
    X_tmp = np.array(X_tmp)
    print(X_tmp.shape)
    Y_total = np.array(Y_total)
    print(Y_total.shape)

    X_tmp = get_pca(input_size, X_tmp)
    X_tmp_size = X_tmp.shape
    print(X_tmp.shape)

    X_verify = X_tmp[0:a[0]]
    X_total = X_tmp[a[0]:X_tmp_size[0]]

    trainX, testX, trainY, testY = model_selection.train_test_split(X_total, Y_total, test_size=0.2, random_state=5)

    trainSize = trainX.shape
    print(testY[1])
    print(Y_total[0:2])
    X_verify = pd.DataFrame(X_verify)
    subclusters = subclusters
    X_verify_total = pd.concat([X_verify, subclusters], axis=1)
    print(X_verify_total.columns)
    print(X_verify_total.shape)
    res_clusters = []
    with tf.Session() as sess:
        init_op = tf.global_variables_initializer()
        sess.run(init_op)
        train_dataset_size = trainSize[0]

        for i in range(training_iterations):
            start = (i * batch_size) % train_dataset_size
            end = min(start + batch_size, train_dataset_size)

            training_loss = sess.run([opt, loss], feed_dict={X: trainX[start:end], Y: trainY[start:end]})
            train_accuracy = accuracy.eval(session=sess, feed_dict={X: trainX[start:end], Y: trainY[start:end]})
            if i % 1000 == 0:
                print("step : %d, training loss = %g " % (i, training_loss[1]))
                print("step : %d, training accuracy = %g " % (i, train_accuracy))

        test_loss = sess.run(loss, feed_dict={X: testX, Y: testY})
        #test_accu = accuracy.eval(session=sess, feed_dict={X: testX, Y: testY})
        print("loss:", test_loss)
        print(accuracy.eval(session=sess, feed_dict={X: testX, Y: testY}))

        preX = testX[1]
        preX = preX.reshape(1, input_size)
        preX = tf.cast(preX, tf.float32)
        print(type(preX))
        print(preX.shape)

        hidden_opt1 = tf.matmul(preX, W1) + B1
        hidden_opt1 = tf.nn.relu(hidden_opt1)
        hidden_opt2 = tf.matmul(hidden_opt1, W2) + B2
        hidden_opt2 = tf.nn.relu(hidden_opt2)
        hidden_opt3 = tf.matmul(hidden_opt2, W3) + B3
        hidden_opt3 = tf.nn.relu(hidden_opt3)
        final_opt = tf.matmul(hidden_opt3, W4) + B4
        a = tf.nn.softmax(final_opt)
        print("result:", sess.run(a))

        subclusters_num = len(set(subclusters.iloc[:,0]))
        print(subclusters_num)
        for num in range(subclusters_num):
            print(num)
            X_verify = X_verify_total[X_verify_total['subclusters'] == str(num)]
            X_verify = X_verify.drop(['subclusters'], axis=1)
            Y_verify = sess.run(y, feed_dict={X: X_verify})
            Y_verify = pd.DataFrame(Y_verify)
            from collections import Counter

            b = Y_verify.shape
            res = []
            for i in range(b[0]):
                a = Y_verify.iloc[i]
                a = a.tolist()
                max_a = max(a)
                idx = a.index(max(a))
                if max_a < 0.35:
                   tmp_res = 'Unknown'
                else:
                   idx = a.index(max(a))
                   tmp_res = dict_idx[idx]
                res.append(tmp_res)
            #pd.DataFrame(res).to_csv(outputpath, index=False)
            print(Counter(res))
            res_clusters.append(max(res, key=res.count))
            #res = pd.DataFrame(res)
    return res_clusters
