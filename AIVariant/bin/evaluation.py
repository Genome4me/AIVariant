###
# Written by Hyeonseong Jeon
# email: jun0605@snu.ac.kr
###

import sys
import os
import pickle 
import tensorflow

#try:
#    import tensorflow
#except ModuleNotFoundError:
#    os.system('pip install tensorflow==2.9.1')

#try:
#    import keras
#except ModuleNotFoundError:
#    os.system('pip install keras==2.9.0')

import tensorflow.keras as keras 
import keras.backend as K 
import math 
import numpy as np 

from keras.models import load_model 

import argparse
import logging
from functools import partial
import pathlib

arch_file = os.path.join('bin', 'architectures', 'AIVariant_model')
checkpoint_file = os.path.join('bin', 'checkpoints', 'aiv_checkpoint')

def get_config():
    parser = argparse.ArgumentParser('config', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument = partial(parser.add_argument, help=' ')

    parser.add_argument('--gpu_id', required=True, help='gpu id')
    parser.add_argument('--test_data_file', required=True, help='file holding test file list and labels')
    parser.add_argument('--prediction_file', required=True, help='prediction file')
    config = parser.parse_args()

    return config

def get_snv_data(files):
    keys = ['tumor_raw_img', 'tumor_processed_img', 'normal_raw_img', 'normal_processed_img', 'tumor_processed_stat_img', 'normal_processed_stat_img', 'lodt', 'lodn', 'unnorm_lodt', 'unnorm_lodn']
    
    X, y = dict(), []
    
    for file in files:
        if '/gp/' in file:
            y.append(1)
        else:
            y.append(0)
            
        with open(file, 'rb') as f:
            data = pickle.load(f)
            for key in keys:
                if key not in X:
                    X[key] = []
                if data[key] is None:
                    X[key].append(0)
                else:
                    X[key].append(data[key])
    
    for key in X:
        X[key] = np.array(X[key])
    y = np.array(y)
    
    return X, y

def get_data(files):
    indel_bases = ['I1', 'I2', 'I3', 'D1', 'D2', 'D3']
    X, y = dict(), []
    
    for file in files:
        
        indel_candidate = any([x in file for x in indel_bases])
        
        if '/gp/' in file:
            if all([x not in os.path.basename(file) for x in indel_bases]):
                y.append([1, 1])
            else:
                alt = [x for x in indel_bases if x in os.path.basename(file)][0]
                b = 2 if alt[0] == 'I' else 3
                s = int(alt[1])
                y.append([b, s])
        else:
            y.append([0, 0])
            
        with open(file, 'rb') as f:
            data = pickle.load(f)
            for key in data:
                if key not in X:
                    X[key] = []
                if data[key] is None:
                    X[key].append(0)
                else:
                    X[key].append(data[key])
    
    for key in X:
        X[key] = np.array(X[key])
    y = np.array(y)
    
    return X, y

def get_snv_label(files):
    y = []
    
    for file in files:
        if '/gp/' in file:
            y.append(1)
        else:
            y.append(0)
    
    y = np.array(y)
    
    return y

def get_label(files):
    indel_bases = ['I1', 'I2', 'I3', 'D1', 'D2', 'D3']
    y = []
    
    for file in files:
        if '/gp/' in file:
            if all([x not in os.path.basename(file) for x in indel_bases]):
                y.append([1, 1])
            else:
                alt = [x for x in indel_bases if x in os.path.basename(file)][0]
                b = 2 if alt[0] == 'I' else 3
                s = int(alt[1])
                y.append([b, s])
        else:
            y.append([0, 0])
    
    y = np.array(y)
    
    return y

def get_candidate_label(files):
    y = []
    
    for file in files:
        alt = os.path.basename(file)[:-4].split(':')[-1]
        if alt in "ACGT":
            y.append([1, 1])
        else:
            if alt[0] == 'I':
                y.append([2, int(alt[1])])
            else:
                y.append([3, int(alt[1])])
                
    return np.array(y)

class AIVSNVDataLoader(keras.utils.Sequence):
    def __init__(self, files, batch_size=32, shuffle=False):
        self.files = np.array(files)
        
        self.batch_size = batch_size
        
        self.shuffle = shuffle
            
        self.on_epoch_end()
            
    def __len__(self):
        return math.ceil(len(self.files)/self.batch_size)
    
    def __getitem__(self, idx):
        indices = self.indices[self.batch_size*idx:self.batch_size*(idx+1)]
        files_batch = self.files[indices]
        return get_snv_data(files_batch)
    
    def on_epoch_end(self):
        self.indices = np.arange(len(self.files))
            
        if self.shuffle:
            np.random.shuffle(self.indices)

class AIVDataLoader(keras.utils.Sequence):
    def __init__(self, files, batch_size=32, shuffle=False):
        self.files = np.array(files)
        
        self.batch_size = batch_size
        
        self.shuffle = shuffle
            
        self.on_epoch_end()
            
    def __len__(self):
        return math.ceil(len(self.files)/self.batch_size)
    
    def __getitem__(self, idx):
        indices = self.indices[self.batch_size*idx:self.batch_size*(idx+1)]
        files_batch = self.files[indices]
        return get_data(files_batch)
    
    def on_epoch_end(self):
        self.indices = np.arange(len(self.files))
            
        if self.shuffle:
            np.random.shuffle(self.indices)

def get_exp_threshold(preds, delta=0.0001):
    hist = np.histogram(preds, bins=np.arange(0, 1+delta, delta))
    
    whole = np.sum(hist[0])
    reference = hist[0][-1]
    
    if hist[0][-1] - hist[0][-2] < 0.02 * whole:
        return 1-delta
    
    THs = hist[1][:-1][::-1]
    vals = hist[0][::-1]
    acc = 0
    for i, th in enumerate(THs):
        if i == 0:
            continue
        acc += vals[i]
        if not acc <= reference or vals[i] - vals[i-1] > 10*delta * whole: 
            return th
        
    return 0

def main():
    config = get_config() 

    # os.environ["CUDA_VISIBLE_DEVICES"] = config.gpu_id 

    model = load_model(arch_file)
    model.load_weights(checkpoint_file)
    # model.compile(loss="sparse_categorical_crossentropy", optimizer=keras.optimizers.Adam())

    # print(config.test_data_file)

    test_snv_files = []
    test_indel_files = [] # We have small indel candidates but do not evaluate these in the current version.
    with open(config.test_data_file) as f:
        # print("file opened")
        indel_bases = ['I1', 'I2', 'I3', 'D1', 'D2', 'D3']
        # print(f.readlines())
        for line in f:
            # print("!!!!!!")
            # print(line)
            # print("!!!!!!")
            filepath = line.strip()
            if not filepath.endswith(".pkl"):
                continue
            if all([x not in os.path.basename(filepath) for x in indel_bases]):
                test_snv_files.append(filepath)
            else:
                test_indel_files.append(filepath)

    # return

    test_dl = AIVSNVDataLoader(test_snv_files) # Only evaluates SNV.

    preds = model.predict(test_dl)

    labels = get_snv_label(test_snv_files)
    # labels = get_label(test_snv_files + test_indel_files)
    # candidate_labels = get_candidate_label(test_snv_files + test_indel_files)

    # preds_snv, preds_indel = preds[:len(test_snv_files)], preds[len(test_snv_files):]
    # labels_snv, labels_indel = labels[:len(test_snv_files)], labels[len(test_snv_files):]
    # candidate_labels_snv, candidate_labels_indel = candidate_labels[:len(test_snv_files)], candidate_labels[len(test_snv_files):]

    threshold = get_exp_threshold(preds)

    out = open(config.prediction_file, "w")

    for i, pred in enumerate(preds):
        if np.argmax(pred, axis=0) == labels[i]:
            score = max(pred) 
            if score < threshold:
                score = 0
            filename = os.path.basename(test_snv_files[i])
            
            chrom, pos, ref, alt = filename.strip()[:-4].split(':')
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\tSNV\t{score}\n")
        else:
            filename = os.path.basename(test_snv_files[i])
            
            chrom, pos, ref, alt = filename.strip()[:-4].split(':')
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\tSNV\t{0}\n")

    # for i, pred in enumerate(preds_indel):
    #     called_label, called_size = np.argmax(pred[:4], axis=0), np.argmax(pred[4:], axis=0)
    #     candidate_label, candidate_size = candidate_labels_indel[i]
        
    #     if called_label > 1 and called_label == candidate_label and called_size == candidate_size:
    #         score = pred[:4][called_label] * pred[4:][called_size]
    #         filename = os.path.basename(test_indel_files[i])
            
    #         chrom, pos, ref, alt = filename.strip()[:-4].split(':')
    #         out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\tINDEL\t{score}\n")

    out.close()

if __name__=="__main__":
    main()
