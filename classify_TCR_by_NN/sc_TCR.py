#!/usr/bin/env python3

# -------------------
# train a model
# -------------------


import os
import sys
import time
import argparse

import random

import numpy as np
import pandas as pd
import tensorflow as tf

sys.path.append(os.getcwd())

from sc_TCR_utils import _utils



parser = argparse.ArgumentParser(
    description=\
    "sc_TCR: a neural network model for cancer-associated single cell TCR prediction",
        prog='sc_TCR')

parser.add_argument('--feature_choice', required=True, type=str,
                    help='what features to use',
                    choices=['alpha_cdr3', 'alpha_cdr3_v', 'beta_cdr3', 'beta_cdr3_v', 'both_cdr3', 'both_cdr3_v'])
parser.add_argument('--model_dir', required=True, type=str,
                    help='the path to the folder to save the trained model')
parser.add_argument('--result_dir', required=True, type=str,
                    help='the path to the folder to save the test results')
parser.add_argument('--score_dir', required=True, type=str,
                    help='the path to the folder to save the test results')
parser.add_argument('--enc_method', type=str, default='one_hot',
                    help='encoding method for amino acid, can be one of one_hot, blosum62, atchley, and pca',
                    choices=['one_hot', 'blosum62', 'atchley', 'pca'])
parser.add_argument('--freq_cutoff', type=int, default=10000,
                    help='any frequency higher than this value will be clipped')
parser.add_argument('--lr', type=float, default=1e-4,
                    help='learning rate, for example, 0.0001')
parser.add_argument('--n_CNN', type=int, default=1,
                    help='number of CNN layers, can be either 1 or 2',
                    choices=[1, 2])
parser.add_argument('--n_dense', type=int, default=2,
                    help='number of dense layers, can be either 1 or 2',
                    choices=[1, 2])
parser.add_argument('--n_units_str', type=str, default='[64,16]',
                    help='a list of sizes of the dense layers, must be in the format of a list '\
                         'with length matching n_dense (with only ",", no space between numbers in the list, '\
                         'if the list has length 2), for example, [64,16]')
parser.add_argument('--dropout_flag', type=str, default='True',
                    help='whether to use dropout or not, can be either True or False',
                    choices=['True', 'False'])
parser.add_argument('--p_dropout', type=float, default=0.2,
                    help='dropout probability, for example, 0.2')
parser.add_argument('--rseed', type=int, default=1000,
                    help='random seed for random')
parser.add_argument('--np_seed', type=int, default=1216,
                    help='random seed for numpy')
parser.add_argument('--tf_seed', type=int, default=2207,
                    help='random seed for tensorflow')



def train(feature_choice, model_dir, result_dir, score_dir, enc_method,
          freq_cutoff, lr, n_CNN, n_dense, n_units_str, dropout_flag, p_dropout,
          rseed, np_seed, tf_seed):

    input_args = locals()
    print("input args are", input_args)

    random.seed(rseed)
    np.random.seed(np_seed)
    tf.random.set_seed(tf_seed)

    patience = 10

    dropout_flag = (dropout_flag == 'True')

    if n_dense > 2:
        sys.exit("Error: number of dense layers>2, which is a case not coded for yet.")

    len_n_units_str = len(n_units_str)
    n_units_str_input = n_units_str[1:(len_n_units_str - 1)].split(',')
    n_units = [int(i) for i in n_units_str_input]
    print("n_units = ", n_units)
    if len(n_units) != n_dense:
        sys.exit("Error: n_dense and n_units do not match.")

    print("feature_choice = ", feature_choice)
    print("model_dir = ", model_dir)
    print("enc_method = ", enc_method)
    print("freq_cutoff = ", freq_cutoff)

    print("lr = ", lr)
    print("n_dense = ", n_dense)
    print("n_units = ", n_units)
    print("dropout_flag = ", dropout_flag)
    print("p_dropout = ", p_dropout)


    setting_name = \
        feature_choice + "_" + enc_method + '_' + str(lr)[2:] + \
        '_CNN' + str(n_CNN) + \
        '_dense' + str(n_dense) + \
        '_n_units_' + '_'.join([str(n) for n in n_units]) + \
        ('_dropout_p_' + str(p_dropout)[2:]) * int(dropout_flag)

    print(setting_name)

    cur_name = feature_choice+"_"+enc_method+"_"+str(n_CNN)+"_CNN"+'_dense'+str(n_dense)+\
               '_n_units_' + '_'.join([str(n) for n in n_units])

    cur_seed_name = setting_name+"_"+"_".join([str(rseed), str(np_seed), str(tf_seed)])

    checkpoint_path = model_dir+"/"+cur_name+"/"+cur_seed_name

    if not os.path.exists(checkpoint_path):
        os.makedirs(checkpoint_path)



    start = time.time()

    ((train_parts, y2_train, w_train),
     (valid_parts, y2_valid),
     (test_parts, y2_test)) = \
        _utils.get_data(enc_method, freq_cutoff)

    print(train_parts[0].shape)
    print(valid_parts[0].shape)
    print(test_parts[0].shape)

    print(train_parts[1].shape)
    print(train_parts[2].shape)
    print(train_parts[3].shape)
    print(train_parts[4].shape)
    print(train_parts[5].shape)
    print(train_parts[6].shape)
    print(train_parts[7].shape)
    print(train_parts[8].shape)
    print(train_parts[9].shape)

    print(y2_train.shape)


    # get the model
    model = _utils.get_model(feature_choice=feature_choice,
                            CDR3_shape=train_parts[0].shape[1:],
                            len_shape=train_parts[2].shape[1:],
                            cdr1_shape=train_parts[4].shape[1:],
                            cdr2_shape=train_parts[6].shape[1:],
                            cdr25_shape=train_parts[8].shape[1:],
                            n_CNN=n_CNN,
                            n_dense=n_dense,
                            n_units=n_units,
                            dropout_flag=dropout_flag,
                            p_dropout=p_dropout)

    model.summary()

    # compile model

    METRICS = [
        #tf.keras.metrics.BinaryAccuracy(name="accuracy"),
        tf.keras.metrics.AUC(name="auc_roc", curve='ROC')
        #tf.keras.metrics.AUC(name="auc_pr", curve='PR')
    ]

    adam_optim = tf.keras.optimizers.Adam(learning_rate=lr)
    model.compile(loss=tf.keras.losses.BinaryCrossentropy(), optimizer=adam_optim,
                  metrics=METRICS)

    # fit model

    sample_weights = np.array(w_train).reshape(len(w_train),)

    callback = tf.keras.callbacks.EarlyStopping(monitor='val_auc_roc',
                                                patience=patience,
                                                mode='max')
    # this check_point saves the model with the best performance on validation data
    check_point = tf.keras.callbacks.ModelCheckpoint(
        checkpoint_path,
        monitor='val_auc_roc',
        verbose=1,
        save_best_only=True,
        mode='max',
        save_weights_only=False
    )

    # validation sample weight seems to be needed for tensorflow version
    # on gpu on server, but not for the version on mac desktop
    # since we are using AUC as the validation metric
    # it doesn't matter whether we use weight or not
    # so it doesn't matter if the weight does not equal the proportion in valid
    y2_valid_list = [y[0] for y in y2_valid]


    model.fit(x=list(train_parts), y=y2_train, sample_weight=sample_weights,
              validation_data=(list(valid_parts), y2_valid),
              callbacks=[callback, check_point],
              epochs=300, batch_size=32)

    end = time.time()
    print("Training the model took around "+str(round((end - start)/60, 1)) + " minutes.")
    print("Finished training.")

    tf.keras.backend.clear_session()

    model = tf.keras.models.load_model(checkpoint_path)


    y2_test_list = [y[0] for y in y2_test]

    results = model.evaluate(test_parts, y2_test)
    yhat_value = model.predict(test_parts)
    yhat_label = (yhat_value>0.5).astype("int32")

    po_ind = [i for i in range(len(y2_test_list)) if y2_test_list[i] == 1]
    yhat_label = yhat_label.reshape(yhat_label.shape[0], 1)
    acc_on_po = yhat_label[po_ind].mean()
    ne_ind = [i for i in range(len(y2_test_list)) if y2_test_list[i] == 0]
    acc_on_ne = 1 - yhat_label[ne_ind].mean()

    df_result = pd.DataFrame(list(zip([results[1]], [acc_on_po], [acc_on_ne])),
                             columns=["auc", "acc_on_po", "acc_on_ne"])

    if not os.path.exists(result_dir+"/"+cur_name):
        os.makedirs(result_dir+"/"+cur_name)

    if not os.path.exists(score_dir+"/"+cur_name):
        os.makedirs(score_dir+"/"+cur_name)

    df_result.to_csv(result_dir+"/"+cur_name+"/"+cur_seed_name +".csv",
                     index=False)

    df_score = pd.DataFrame(yhat_value, columns = ['score'])
    df_score.to_csv(score_dir+"/"+cur_name+"/"+cur_seed_name +".csv",
                     index=False)



if __name__ == "__main__":
    args = parser.parse_args()
    print(args)
    train(**vars(args))
