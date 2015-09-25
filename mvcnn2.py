# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# add to kfkd.py
from lasagne import layers
from lasagne.updates import nesterov_momentum
from nolearn.lasagne import NeuralNet
import numpy as np
import theano.tensor as T
from nolearn.lasagne import BatchIterator
from theano.sandbox.neighbours import neibs2images

from lasagne.nonlinearities import tanh
from lasagne.nonlinearities import rectify
import pickle
import sys
from sklearn.metrics import mean_squared_error as mse
from sklearn.metrics import precision_score
import os
import urllib
import gzip
import cPickle
from IPython.display import Image as IPImage
from PIL import Image


# read in data as a pickle file (example is shown below for mnist)

# fname = 'mnist/mnist.pkl.gz'
# if not os.path.isfile(fname):
#     testfile = urllib.URLopener()
#     testfile.retrieve("http://deeplearning.net/data/mnist/mnist.pkl.gz", fname)
# f = gzip.open(fname, 'rb')
# train_set, valid_set, test_set = cPickle.load(f)
# f.close()
# X, y = train_set
# X = np.rint(X * 256).astype(np.int).reshape((-1, 1, 28, 28))  # convert to (0,255) int range (we'll do our own scaling)
# mu, sigma = np.mean(X.flatten()), np.std(X.flatten())

# X_train = X.astype(np.float64)
# X_train = (X_train - mu) / sigma
# X_train = X_train.astype(np.float32)

# # we need our target to be 1 dimensional
# X_out = X_train.reshape((X_train.shape[0], -1))


conv_filters = 32
filter_sizes = 7
epochs = 20
encode_size = 40
num_classes=10

# 3 multi-view inputs. These go in parallel to the N (e.g. N=3) CNNs shown below
l_in1 = layers.InputLayer((None, 1, 28, 28))
l_in2 = layers.InputLayer((None, 1, 28, 28))
l_in3 = layers.InputLayer((None, 1, 28, 28))

# parallel CNNs for each input (each CNN has two convolutional layers)
l_conv1 = layers.Conv2DLayer(l_in1, num_filters=conv_filters, filter_size = (filter_sizes, filter_sizes), pad="valid", nonlinearity=rectify, W=init.GlorotUniform(), b=init.Constant(0.))
l_conv2 = layers.Conv2DLayer(l_in2, num_filters=conv_filters, filter_size = (filter_sizes, filter_sizes), pad="valid", nonlinearity=rectify, W=init.GlorotUniform(), b=init.Constant(0.))
l_conv3 = layers.Conv2DLayer(l_in3, num_filters=conv_filters, filter_size = (filter_sizes, filter_sizes), pad="valid", nonlinearity=rectify, W=init.GlorotUniform(), b=init.Constant(0.))

l_pool1 = layers.MaxPool2DLayer(l_conv1 pool_size=(2,2))
l_pool2 = layers.MaxPool2DLayer(l_conv2, pool_size=(2,2))
l_pool3 = layers.MaxPool2DLayer(l_conv3, pool_size=(2,2))

l_conv21 = layers.Conv2DLayer(l_pool1, num_filters=conv_filters, filter_size = (filter_sizes, filter_sizes), pad="valid", nonlinearity=rectify, W=init.GlorotUniform(), b=init.Constant(0.))
l_conv22 = layers.Conv2DLayer(l_pool2, num_filters=conv_filters, filter_size = (filter_sizes, filter_sizes), pad="valid", nonlinearity=rectify, W=init.GlorotUniform(), b=init.Constant(0.))
l_conv23 = layers.Conv2DLayer(l_pool3, num_filters=conv_filters, filter_size = (filter_sizes, filter_sizes), pad="valid", nonlinearity=rectify, W=init.GlorotUniform(), b=init.Constant(0.))

l_pool21 = layers.MaxPool2DLayer(l_conv21 pool_size=(2,2))
l_pool22 = layers.MaxPool2DLayer(l_conv22, pool_size=(2,2))
l_pool23 = layers.MaxPool2DLayer(l_conv23, pool_size=(2,2))

# common layer concatenating all the parallel N CNNs
l_comm = lasagne.layers.ConcatLayer([l_pool21, l_pool22, l_pool23])

# a common CNN now feeds forward to the output layer
l_comm_conv1 = layers.Conv2DLayer(l_comm, num_filters=conv_filters, filter_size = (filter_sizes, filter_sizes), pad="valid", nonlinearity=None, W=init.GlorotUniform(), b=init.Constant(0.))
l_comm_pool1 = layers.MaxPool2DLayer(l_conv1 pool_size=(2,2))

l_comm_conv2 = layers.Conv2DLayer(l_comm_pool1, num_filters=conv_filters, filter_size = (filter_sizes, filter_sizes), pad="valid", nonlinearity=None, W=init.GlorotUniform(), b=init.Constant(0.))
l_comm_pool2 = layers.MaxPool2DLayer(l_comm_conv2 pool_size=(2,2))

# output softmax layer for classification into num_classes classes
l_out = lasagne.layers.DenseLayer(
        l_comm_pool2, num_units=num_classes,
        nonlinearity=lasagne.nonlinearities.softmax)