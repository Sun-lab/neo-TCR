

import numpy as np
import pandas as pd

import random

import pkg_resources
#from importlib import resources

from collections import defaultdict

from sklearn.preprocessing import OneHotEncoder

from tensorflow.keras.activations import relu
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense, Flatten, Conv1D, MaxPooling1D
from tensorflow.keras.layers import Reshape, Dropout, concatenate

from tensorflow.keras.callbacks import Callback




def prepare_encoders(enc_method):



    # hard encode sequence length limit
    unique_TCR_aa_lens = list(range(1, 28))
    CDR3len_enc_template = \
        np.array(unique_TCR_aa_lens).reshape(len(unique_TCR_aa_lens), 1)
    CDR3len_enc = OneHotEncoder().fit(CDR3len_enc_template)

    # build encoder for amino acids in CDR3
    # one_hot_encoding first
    if enc_method == "one_hot":
        AA_SYMOLS = ['A', 'R', 'N', 'D', 'C',
                     'Q', 'E', 'G', 'H', 'I',
                     'L', 'K', 'M', 'F', 'P',
                     'S', 'T', 'W', 'Y', 'V']
        AA_SYMOLS.sort()
        # now move on to encode the v gene and CDR3 part
        # add '.' to the aa name list for the purpose of padding short CDR3 sequences
        # encoder for CDR3
        # need to use max length of CDR3 aa seq in current file
        AA_SYMOLS_pad = AA_SYMOLS + ['.']
        CDR3_enc_template = np.array([[aa] * 27 for aa in AA_SYMOLS_pad])
        CDR3_enc = OneHotEncoder().fit(CDR3_enc_template)
        # encoders for cdr1, cdr2, cdr2.5
        # due to the fact that the cdr1, cdr2 sequence contains '*', the CDR3_enc below is
        # adjusted to add a '*' option to make it work
        # for cdr1, cdr2, cdr25 as weill
        AA_SYMOLS_pad_star = AA_SYMOLS + ['.'] + ['*']
        cdr1_template = np.array([[aa] * 12 for aa in AA_SYMOLS_pad_star])
        cdr1_enc = OneHotEncoder().fit(cdr1_template)
        cdr2_template = np.array([[aa] * 10 for aa in AA_SYMOLS_pad_star])
        cdr2_enc = OneHotEncoder().fit(cdr2_template)
        cdr25_template = np.array([[aa] * 6 for aa in AA_SYMOLS_pad])
        cdr25_enc = OneHotEncoder().fit(cdr25_template)
    # or blosum62
    elif enc_method == "blosum62":
        # define encoder for CDR3
        # this encoder also works for cdr25
        # it does not work for cdr1, cdr2 due to that some sequence in cdr1, cdr2 contains
        # a character "*"
        def cdr_enc(cdr_list):

            blosum62_file = 'data/for_encoders/blosum62.csv'
            blosum62_stream = pkg_resources.resource_stream(__name__, blosum62_file)
            #with resources.path("DePTH.data", blosum62_file) as blosum62_stream:
            blosum62_matrix = pd.read_csv(blosum62_stream, sep=',', header=0)

            blosum62_array = np.array(blosum62_matrix)
            blosum62_dict = defaultdict(list)
            for i in range(20):
                blosum62_dict[blosum62_array[i][0]] = \
                    blosum62_array[i][1:].tolist() + [-4]
            blosum62_dict['.'] = [-4 for _ in range(20)] + [1]
            return [[blosum62_dict[aa] for aa in cdr] for cdr in cdr_list]

        # define an encoder specifically for cdr1 and cdr2
        # since some sequence in cdr1, cdr2 contains a character "*"
        def cdr2_enc(cdr2_list):

            blosum62_file = 'data/for_encoders/blosum62.csv'
            blosum62_stream = pkg_resources.resource_stream(__name__, blosum62_file)
            #with resources.path("DePTH.data", blosum62_file) as blosum62_stream:
            blosum62_matrix = pd.read_csv(blosum62_stream, sep=',', header=0)

            blosum62_array = np.array(blosum62_matrix)
            blosum62_dict = defaultdict(list)
            for i in range(20):
                blosum62_dict[blosum62_array[i][0]] = \
                    blosum62_array[i][1:].tolist() + [-4]
            blosum62_dict['.'] = [-4 for _ in range(20)] + [1]
            blosum62_dict['*'] = [-4 for _ in range(20)] + [1]
            return [[blosum62_dict[aa] for aa in cdr2] for cdr2 in cdr2_list]
    # or atchley
    elif enc_method == "atchley":
        # define encoder for CDR3
        # this encoder also works for cdr25
        # it does not work for cdr1, cdr2 due to that some sequence in cdr1, cdr2 contains
        # a charater "*"
        def cdr_enc(cdr_list):

            atchley_file = 'data/for_encoders/Atchley_factors.csv'
            atchley_stream = pkg_resources.resource_stream(__name__, atchley_file)
            #with resources.path("DePTH.data", atchley_file) as atchley_stream:
            atchley_matrix = pd.read_csv(atchley_stream, sep=',', header=0)

            atchley_array = np.array(atchley_matrix)
            atchley_dict = defaultdict(list)
            for i in range(20):
                atchley_dict[atchley_array[i][0]] = atchley_array[i][1:].tolist()
            atchley_dict['.'] = [0 for _ in range(5)]
            return [[atchley_dict[aa] for aa in cdr] for cdr in cdr_list]

        # define an encoder specifically for cdr1 and cdr2
        # since some sequence in cdr1, cdr2 contains a character "*"
        def cdr2_enc(cdr2_list):

            atchley_file = 'data/for_encoders/Atchley_factors.csv'
            atchley_stream = pkg_resources.resource_stream(__name__, atchley_file)
            #with resources.path("DePTH.data", atchley_file) as atchley_stream:
            atchley_matrix = pd.read_csv(atchley_stream, sep=',', header=0)

            atchley_array = np.array(atchley_matrix)
            atchley_dict = defaultdict(list)
            for i in range(20):
                atchley_dict[atchley_array[i][0]] = atchley_array[i][1:].tolist()
            atchley_dict['.'] = [0 for _ in range(5)]
            atchley_dict['*'] = [0 for _ in range(5)]
            return [[atchley_dict[aa] for aa in cdr2] for cdr2 in cdr2_list]
    # AAidx_PCA
    else:
        # define encoder for CDR3
        # this encoder also works for cdr25
        # it does not work for cdr1, cdr2 due to that some sequence in cdr1, cdr2 contains
        # a character "*"
        def cdr_enc(cdr_list):

            pca_file = 'data/for_encoders/AAidx_PCA.csv'
            pca_stream = pkg_resources.resource_stream(__name__, pca_file)
            #with resources.path("DePTH.data", pca_file) as pca_stream:
            pca_matrix = pd.read_csv(pca_stream, sep=',', header=0)

            pca_array = np.array(pca_matrix)
            pca_dict = defaultdict(list)
            for i in range(20):
                pca_dict[pca_array[i][0]] = pca_array[i][1:].tolist()
            pca_dict['.'] = [0 for _ in range(15)]
            return [[pca_dict[aa] for aa in cdr] for cdr in cdr_list]

        # define an encoder specifically for cdr1, cdr2
        # since some sequence in cdr1, cdr2 contains a charater "*"
        def cdr2_enc(cdr2_list):

            pca_file = 'data/for_encoders/AAidx_PCA.csv'
            pca_stream = pkg_resources.resource_stream(__name__, pca_file)
            #with resources.path("DePTH.data", pca_file) as pca_stream:
            pca_matrix = pd.read_csv(pca_stream, sep=',', header=0)

            pca_array = np.array(pca_matrix)
            pca_dict = defaultdict(list)
            for i in range(20):
                pca_dict[pca_array[i][0]] = pca_array[i][1:].tolist()
            pca_dict['.'] = [0 for _ in range(15)]
            pca_dict['*'] = [0 for _ in range(15)]
            return [[pca_dict[aa] for aa in cdr2] for cdr2 in cdr2_list]


    if enc_method == "one_hot":
        return (CDR3len_enc, CDR3_enc, cdr1_enc, cdr2_enc, cdr25_enc)
    else:
        return (CDR3len_enc, cdr_enc, cdr2_enc, cdr2_enc, cdr_enc)


# define padding function for later use
def pad_cdr3(cdr3, max_length, pad_letter):
    if len(cdr3) == max_length:
        return cdr3
    else:
        gap_start = min(6, 3 + (len(cdr3) - 5) // 2)
        return cdr3[:gap_start] + [pad_letter for _ in range(max_length - len(cdr3))] + cdr3[gap_start:]


def encode(items_ori, enc_method, CDR3len_enc, CDR3_enc,
            cdr1_enc, cdr2_enc, cdr25_enc):

    # load info for cdr1, cdr2, cdr25 encoding
    V_info_file = 'data/for_encoders/combo_xcr.tsv'
    V_stream = pkg_resources.resource_stream(__name__, V_info_file)
    #with resources.path("DePTH.data", V_info_file) as V_stream:
    V_info = pd.read_csv(V_stream, sep='\t')

    V_sub_info = V_info.loc[(V_info.organism == 'human')
                            & (V_info.chain.isin(['A', 'B']))
                            & (V_info.region == 'V')]
    cdr1_dict = defaultdict(str)
    cdr2_dict = defaultdict(str)
    cdr25_dict = defaultdict(str)

    for allele, cdrs in zip(V_sub_info.id.tolist(), V_sub_info.cdrs.tolist()):
        cdr_seqs = cdrs.split(";")
        cdr1_dict[allele] = cdr_seqs[0]
        cdr2_dict[allele] = cdr_seqs[1]
        cdr25_dict[allele] = cdr_seqs[2]

    # define a function to transform va, vb
    special_dict = defaultdict(str)

    special_dict['TRAV14DV4'] = 'TRAV14/DV4*01'
    special_dict['TRAV23DV6'] = 'TRAV23/DV6*01'
    special_dict['TRAV29DV5'] = 'TRAV29/DV5*01'
    special_dict['TRAV36DV7'] = 'TRAV36/DV7*01'
    special_dict['TRAV38-2DV8'] = 'TRAV38-2/DV8*01'

    def trans_vab(raw_v):
        if raw_v not in special_dict:
            return raw_v+"*01"
        else:
            return special_dict[raw_v]

    # encode
    # CDR3

    CDR3_a_list = [pad_cdr3(list(tcr.split(",")[2]), 27, '.')
                   for tcr in items_ori]
    CDR3_b_list = [pad_cdr3(list(tcr.split(",")[3]), 27, '.')
                   for tcr in items_ori]
    # cdr1, cdr2, cdr25
    V_allele_trans_a = [trans_vab(tcr.split(",")[0]) for tcr in items_ori]
    V_allele_trans_b = [trans_vab(tcr.split(",")[1]) for tcr in items_ori]

    cdr1_a_list = [list(cdr1_dict[trans]) for trans in V_allele_trans_a]
    cdr2_a_list = [list(cdr2_dict[trans]) for trans in V_allele_trans_a]
    cdr25_a_list = [list(cdr25_dict[trans]) for trans in V_allele_trans_a]

    cdr1_b_list = [list(cdr1_dict[trans]) for trans in V_allele_trans_b]
    cdr2_b_list = [list(cdr2_dict[trans]) for trans in V_allele_trans_b]
    cdr25_b_list = [list(cdr25_dict[trans]) for trans in V_allele_trans_b]

    # length of cdr3
    CDR3_len_a_list = [[len(tcr.split(",")[2])] for tcr in items_ori]
    CDR3_len_a = np.array(CDR3_len_a_list)
    CDR3len_a_encoded = CDR3len_enc.transform(CDR3_len_a).toarray()

    CDR3_len_b_list = [[len(tcr.split(",")[3])] for tcr in items_ori]
    CDR3_len_b = np.array(CDR3_len_b_list)
    CDR3len_b_encoded = CDR3len_enc.transform(CDR3_len_b).toarray()

    # encoding for amino acis in CDR3
    if enc_method == "one_hot":
        # CDR3
        CDR3_a_array = np.array(CDR3_a_list)
        CDR3_a_encoded = CDR3_enc.transform(CDR3_a_array).toarray()
        CDR3_a_encoded = CDR3_a_encoded.reshape(len(items_ori), 27, 21)
        CDR3_b_array = np.array(CDR3_b_list)
        CDR3_b_encoded = CDR3_enc.transform(CDR3_b_array).toarray()
        CDR3_b_encoded = CDR3_b_encoded.reshape(len(items_ori), 27, 21)
        # cdr1, cdr2, cdr25
        cdr1_a_array = np.array(cdr1_a_list)
        cdr1_a_encoded = cdr1_enc.transform(cdr1_a_array).toarray()
        cdr1_a_encoded = cdr1_a_encoded.reshape(len(items_ori), 12, 22)
        cdr1_b_array = np.array(cdr1_b_list)
        cdr1_b_encoded = cdr1_enc.transform(cdr1_b_array).toarray()
        cdr1_b_encoded = cdr1_b_encoded.reshape(len(items_ori), 12, 22)

        cdr2_a_array = np.array(cdr2_a_list)
        cdr2_a_encoded = cdr2_enc.transform(cdr2_a_array).toarray()
        cdr2_a_encoded = cdr2_a_encoded.reshape(len(items_ori), 10, 22)
        cdr2_b_array = np.array(cdr2_b_list)
        cdr2_b_encoded = cdr2_enc.transform(cdr2_b_array).toarray()
        cdr2_b_encoded = cdr2_b_encoded.reshape(len(items_ori), 10, 22)

        cdr25_a_array = np.array(cdr25_a_list)
        cdr25_a_encoded = cdr25_enc.transform(cdr25_a_array).toarray()
        cdr25_a_encoded = cdr25_a_encoded.reshape(len(items_ori), 6, 21)
        cdr25_b_array = np.array(cdr25_b_list)
        cdr25_b_encoded = cdr25_enc.transform(cdr25_b_array).toarray()
        cdr25_b_encoded = cdr25_b_encoded.reshape(len(items_ori), 6, 21)
    # blosum62, atchley, pca
    else:
        CDR3_a_encoded = np.array(CDR3_enc(CDR3_a_list))
        CDR3_b_encoded = np.array(CDR3_enc(CDR3_b_list))
        cdr1_a_encoded = np.array(cdr1_enc(cdr1_a_list))
        cdr1_b_encoded = np.array(cdr1_enc(cdr1_b_list))
        cdr2_a_encoded = np.array(cdr2_enc(cdr2_a_list))
        cdr2_b_encoded = np.array(cdr2_enc(cdr2_b_list))
        cdr25_a_encoded = np.array(cdr25_enc(cdr25_a_list))
        cdr25_b_encoded = np.array(cdr25_enc(cdr25_b_list))

    # return encoded data
    return CDR3_a_encoded, CDR3_b_encoded, \
           CDR3len_a_encoded, CDR3len_b_encoded, \
           cdr1_a_encoded, cdr1_b_encoded, \
           cdr2_a_encoded, cdr2_b_encoded, \
           cdr25_a_encoded, cdr25_b_encoded


def get_data(enc_method, freq_cutoff):

    folder_path = "./data_for_model/"

    df_pos_train = pd.read_csv(folder_path+"train_pos.csv", header=0)
    df_pos_valid = pd.read_csv(folder_path+"valid_pos.csv", header=0)
    df_pos_test = pd.read_csv(folder_path+"test_pos.csv", header=0)

    df_neg_train = pd.read_csv(folder_path+"train_neg.csv", header=0)
    df_neg_valid = pd.read_csv(folder_path+"valid_neg.csv", header=0)
    df_neg_test = pd.read_csv(folder_path+"test_neg.csv", header=0)

    train_pos_f = [tuple(row) for row in df_pos_train.to_numpy()]
    train_neg_f = [tuple(row) for row in df_neg_train.to_numpy()]

    valid_pos_f = [tuple(row) for row in df_pos_valid.to_numpy()]
    valid_neg_f = [tuple(row) for row in df_neg_valid.to_numpy()]

    test_pos_f = [tuple(row) for row in df_pos_test.to_numpy()]
    test_neg_f = [tuple(row) for row in df_neg_test.to_numpy()]

    train_pos = [x for (x, f) in train_pos_f]
    train_neg = [x for (x, f) in train_neg_f]

    valid_pos = [x for (x, f) in valid_pos_f]
    valid_neg = [x for (x, f) in valid_neg_f]

    test_pos = [x for (x, f) in test_pos_f]
    test_neg = [x for (x, f) in test_neg_f]

    temp_pos_w_train = [f if f < freq_cutoff else freq_cutoff for (x, f) in train_pos_f]
    temp_neg_w_train = [f if f < freq_cutoff else freq_cutoff for (x, f) in train_neg_f]

    temp_pos_w_valid = [f if f < freq_cutoff else freq_cutoff for (x, f) in valid_pos_f]
    temp_neg_w_valid = [f if f < freq_cutoff else freq_cutoff for (x, f) in valid_neg_f]

    temp_pos_w_test = [f if f < freq_cutoff else freq_cutoff for (x, f) in test_pos_f]
    temp_neg_w_test = [f if f < freq_cutoff else freq_cutoff for (x, f) in test_neg_f]

    sum_pos_w_train = sum(temp_pos_w_train)
    sum_neg_w_train = sum(temp_neg_w_train)

    sum_pos_w_valid = sum(temp_pos_w_valid)
    sum_neg_w_valid = sum(temp_neg_w_valid)

    sum_pos_w_test = sum(temp_pos_w_test)
    sum_neg_w_test = sum(temp_neg_w_test)

    w_ratio_train = sum_neg_w_train/sum_pos_w_train
    w_ratio_valid = sum_neg_w_valid/sum_pos_w_valid
    w_ratio_test = sum_neg_w_test/sum_pos_w_test

    train_pos_w = [f*w_ratio_train for f in temp_pos_w_train]
    train_neg_w = temp_neg_w_train

    valid_pos_w = [f*w_ratio_valid for f in temp_pos_w_valid]
    valid_neg_w = temp_neg_w_valid

    test_pos_w = [f*w_ratio_test for f in temp_pos_w_test]
    test_neg_w = temp_neg_w_test

    # combine pos and neg tcrs together, for training/validation/test separately
    train_combined = train_pos + train_neg
    valid_combined = valid_pos + valid_neg
    test_combined = test_pos + test_neg

    # create labels corresponding to the pos and neg pairs
    train_y = np.array([1 for _ in range(len(train_pos))] +
                       [0 for _ in range(len(train_neg))]).reshape(len(train_pos) + len(train_neg), 1)
    valid_y = np.array([1 for _ in range(len(valid_pos))] +
                       [0 for _ in range(len(valid_neg))]).reshape(len(valid_pos) + len(valid_neg), 1)
    test_y = np.array([1 for _ in range(len(test_pos))] +
                      [0 for _ in range(len(test_neg))]).reshape(len(test_pos) + len(test_neg), 1)

    # combined weights
    train_w_raw = train_pos_w + train_neg_w
    valid_w_raw = valid_pos_w + valid_neg_w
    test_w_raw = test_pos_w + test_neg_w

    # shuffle pairs for training/validation/test separately
    # use the indexes for shuffling to make the pair and label consistent
    # no shuffle on the test data
    train_ind_list = list(range(len(train_combined)))
    valid_ind_list = list(range(len(valid_combined)))
    test_ind_list = list(range(len(test_combined)))
    random.shuffle(train_ind_list)
    random.shuffle(valid_ind_list)
    #random.shuffle(test_ind_list)
    # get shuffled TCRs
    # no shuffle on the test data
    train_ori = [train_combined[i] for i in train_ind_list]
    valid_ori = [valid_combined[i] for i in valid_ind_list]
    test_ori = [test_combined[i] for i in test_ind_list]
    # get shuffled labels
    # no shuffle on the test data
    train_y = train_y[train_ind_list]
    valid_y = valid_y[valid_ind_list]
    test_y = test_y[test_ind_list]
    # get shuffled weights
    # no shuffle on the test data
    train_w = [train_w_raw[i] for i in train_ind_list]
    valid_w = [valid_w_raw[i] for i in valid_ind_list]
    test_w = [test_w_raw[i] for i in test_ind_list]
    # get the elements for encoding
    (CDR3len_enc, CDR3_enc, cdr1_enc, cdr2_enc, cdr25_enc) = prepare_encoders(enc_method)

    print("first two training pairs: ")
    print(train_ori[:2])
    print("first two validation pairs: ")
    print(valid_ori[:2])
    print("first two test pairs: ")
    print(test_ori[:2])

    # get encoded pairs
    components_train = encode(train_ori, enc_method, CDR3len_enc, CDR3_enc,
                              cdr1_enc, cdr2_enc, cdr25_enc)
    components_valid = encode(valid_ori, enc_method, CDR3len_enc, CDR3_enc,
                              cdr1_enc, cdr2_enc, cdr25_enc)
    components_test = encode(test_ori, enc_method, CDR3len_enc, CDR3_enc,
                              cdr1_enc, cdr2_enc, cdr25_enc)

    return ((components_train, train_y, train_w),
            (components_valid, valid_y),
            (components_test, test_y))



def get_model(feature_choice, CDR3_shape, len_shape,
              cdr1_shape, cdr2_shape, cdr25_shape,
              n_CNN, n_dense, n_units, dropout_flag, p_dropout):

    # Define input layers
    CDR3_a_input = Input(CDR3_shape)
    CDR3_b_input = Input(CDR3_shape)

    len_a_input = Input(len_shape)
    len_b_input = Input(len_shape)

    cdr1_a_input = Input(cdr1_shape)
    cdr1_b_input = Input(cdr1_shape)

    cdr2_a_input = Input(cdr2_shape)
    cdr2_b_input = Input(cdr2_shape)

    cdr25_a_input = Input(cdr25_shape)
    cdr25_b_input = Input(cdr25_shape)

    cdr1_a_reshape = Reshape((cdr1_shape[0] * cdr1_shape[1],),
                              input_shape=cdr1_shape)(cdr1_a_input)
    cdr1_b_reshape = Reshape((cdr1_shape[0] * cdr1_shape[1],),
                              input_shape=cdr1_shape)(cdr1_b_input)

    cdr2_a_reshape = Reshape((cdr2_shape[0] * cdr2_shape[1],),
                              input_shape=cdr2_shape)(cdr2_a_input)
    cdr2_b_reshape = Reshape((cdr2_shape[0] * cdr2_shape[1],),
                              input_shape=cdr2_shape)(cdr2_b_input)

    cdr25_a_reshape = Reshape((cdr25_shape[0] * cdr25_shape[1],),
                               input_shape=cdr25_shape)(cdr25_a_input)
    cdr25_b_reshape = Reshape((cdr25_shape[0] * cdr25_shape[1],),
                               input_shape=cdr25_shape)(cdr25_b_input)


    if feature_choice in ["beta_cdr3", "beta_cdr3_v"]:
        CDR3_input = CDR3_b_input
        # construct CDR3_branches
        CDR3_branch = Conv1D(filters=8, kernel_size=2, activation=relu,
                             input_shape=CDR3_shape, name='Conv_CDR3_1')(CDR3_input)
        CDR3_branch = MaxPooling1D(pool_size=2, strides=1, padding='valid',
                                   name='MaxPooling_CDR3_1')(CDR3_branch)
        if n_CNN==2:
            CDR3_branch = Conv1D(filters=16, kernel_size=2, activation=relu, \
                                 input_shape=(25, 8), name='Conv_CDR3_2')(CDR3_branch)
            CDR3_branch = MaxPooling1D(pool_size=2, strides=1, padding='valid', \
                                       name='MaxPooling_CDR3_2')(CDR3_branch)

        CDR3_flatten = Flatten(name='Flatten_CDR3')(CDR3_branch)

        CDR3_reshape = Reshape((CDR3_shape[0] * CDR3_shape[1],),
                               input_shape=CDR3_shape)(CDR3_input)

        CDR3_inter_layer = concatenate([CDR3_flatten, CDR3_reshape], axis=-1)

    elif feature_choice in ["alpha_cdr3", "alpha_cdr3_v"]:
        CDR3_input = CDR3_a_input
        # construct CDR3_branches
        CDR3_branch = Conv1D(filters=8, kernel_size=2, activation=relu,
                             input_shape=CDR3_shape, name='Conv_CDR3_1')(CDR3_input)
        CDR3_branch = MaxPooling1D(pool_size=2, strides=1, padding='valid',
                                   name='MaxPooling_CDR3_1')(CDR3_branch)
        if n_CNN==2:
            CDR3_branch = Conv1D(filters=16, kernel_size=2, activation=relu, \
                                 input_shape=(25, 8), name='Conv_CDR3_2')(CDR3_branch)
            CDR3_branch = MaxPooling1D(pool_size=2, strides=1, padding='valid', \
                                       name='MaxPooling_CDR3_2')(CDR3_branch)

        CDR3_flatten = Flatten(name='Flatten_CDR3')(CDR3_branch)

        CDR3_reshape = Reshape((CDR3_shape[0] * CDR3_shape[1],),
                               input_shape=CDR3_shape)(CDR3_input)

        CDR3_inter_layer = concatenate([CDR3_flatten, CDR3_reshape], axis=-1)

    elif feature_choice in ["both_cdr3", "both_cdr3_v"]:
        CDR3_input = concatenate([CDR3_a_input, CDR3_b_input], axis=1)
        CDR3_combined_shape = (2*CDR3_shape[0], CDR3_shape[1])
        # construct CDR3_branches
        CDR3_branch = Conv1D(filters=16, kernel_size=2, activation=relu,
                             input_shape=CDR3_combined_shape, name='Conv_CDR3_1')(CDR3_input)
        CDR3_branch = MaxPooling1D(pool_size=2, strides=1, padding='valid',
                                   name='MaxPooling_CDR3_1')(CDR3_branch)
        if n_CNN==2:
            CDR3_branch = Conv1D(filters=32, kernel_size=2, activation=relu, \
                                 input_shape=(52, 16), name='Conv_CDR3_2')(CDR3_branch)
            CDR3_branch = MaxPooling1D(pool_size=2, strides=1, padding='valid', \
                                       name='MaxPooling_CDR3_2')(CDR3_branch)

        CDR3_flatten = Flatten(name='Flatten_CDR3')(CDR3_branch)
        CDR3_a_reshape = Reshape((CDR3_shape[0] * CDR3_shape[1],),
                                  input_shape=CDR3_shape)(CDR3_a_input)
        CDR3_b_reshape = Reshape((CDR3_shape[0] * CDR3_shape[1],),
                                  input_shape=CDR3_shape)(CDR3_b_input)
        CDR3_inter_layer = concatenate([CDR3_flatten, CDR3_a_reshape, CDR3_b_reshape], axis=-1)

    else:
        sys.exit("The value of feature_choice is not one of the allowed ones.")


    # concatenate parts together

    if feature_choice == "beta_cdr3":

        TCR_combined = concatenate([len_b_input, CDR3_inter_layer])
        #inter_layer = Dense(32, activation=relu)(TCR_combined)

    elif feature_choice == "alpha_cdr3":

        TCR_combined = concatenate([len_a_input, CDR3_inter_layer])
        #inter_layer = Dense(32, activation=relu)(TCR_combined)

    elif feature_choice == "beta_cdr3_v":

        TCR_combined = concatenate([len_b_input, CDR3_inter_layer,
                                    cdr1_b_reshape, cdr2_b_reshape, cdr25_b_reshape])
        #inter_layer = Dense(64, activation=relu)(TCR_combined)

    elif feature_choice == "alpha_cdr3_v":

        TCR_combined = concatenate([len_a_input, CDR3_inter_layer,
                                    cdr1_a_reshape, cdr2_a_reshape, cdr25_a_reshape])
        #inter_layer = Dense(64, activation=relu)(TCR_combined)

    elif feature_choice == "both_cdr3":

        TCR_combined = concatenate([len_a_input, len_b_input, CDR3_inter_layer])
        #inter_layer = Dense(128, activation=relu)(TCR_combined)

    elif feature_choice == "both_cdr3_v":

        TCR_combined = concatenate([len_a_input, len_b_input, CDR3_inter_layer,
                                    cdr1_a_reshape, cdr1_b_reshape,
                                    cdr2_a_reshape, cdr2_b_reshape,
                                    cdr25_a_reshape, cdr25_b_reshape])
        #inter_layer = Dense(128, activation=relu)(TCR_combined)
    else:
        sys.exit("The value of feature_choice is not one of the allowed ones.")

    # move on to see how many dense layers we want
    # and whether we want a dropout layer
    if n_dense == 1:
        if not dropout_flag:
            last_layer = Dense(n_units[0], activation=relu)(TCR_combined)
        else:
            dense_layer = Dense(n_units[0], activation=relu)(TCR_combined)
            last_layer = Dropout(p_dropout)(dense_layer)
    else:
        if not dropout_flag:
            first_dense = Dense(n_units[0], activation=relu)(TCR_combined)
            last_layer = Dense(n_units[1], activation=relu)(first_dense)
        else:
            first_dense = Dense(n_units[0], activation=relu)(TCR_combined)
            dropout_layer = Dropout(p_dropout)(first_dense)
            last_layer = Dense(n_units[1], activation=relu)(dropout_layer)
    # final output layer
    output = Dense(1, activation='sigmoid', name='output')(last_layer)
    # build the model
    model = Model(inputs=[CDR3_a_input, CDR3_b_input,
                          len_a_input, len_b_input,
                          cdr1_a_input, cdr1_b_input,
                          cdr2_a_input, cdr2_b_input,
                          cdr25_a_input, cdr25_b_input], outputs=output)
    return model
