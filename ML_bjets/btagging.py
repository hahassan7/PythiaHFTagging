#!/usr/bin/env python
from sklearn.model_selection import train_test_split
import btagging_helpers, btagging_models
import pandas as pd
import sys, os
import tensorflow as tf

#### Settings
recreateData            = True
checkData               = False
trainModels             = True
evalModel               = True
predictData             = False
testSize                = 0.2
trainSize               = 0.8
btagging_helpers.gTarget = 2 # 2 = bjets, 1 =cjets

# tf.debugging.set_log_device_placement(True)

data_MC = pd.DataFrame()
dataname = sys.argv[1]
nJetsSample = int(sys.argv[2])


try:
    os.makedirs("Data/{}".format(dataname), exist_ok=True)
    os.makedirs("Results/{}".format(dataname), exist_ok=True)
    os.makedirs("Models/{}".format(dataname), exist_ok=True)
except OSError as error:
    print(f"Error creating directory : {error}")

####### Create filtered dataset, create features etc.
if trainModels or evalModel:
  data_MC = btagging_helpers.LoadInputData_MC(recreate=recreateData, dataname=dataname, nJets=nJetsSample)
  print("Finished loading MC")

#  data_LHC18b8 = btagging_helpers.LoadInputData_LHC18b8(recreate=recreateData)

if checkData:
  print("CheckData requested")
  btagging_helpers.CheckData(data_MC, dataname)
  print(data_MC.corr()['jetFlavor'].sort_values(ascending=False))

####### Split dataset into training and validation
if trainModels or evalModel:
  toy_train, toy_test = train_test_split(data_MC, test_size=testSize, train_size =trainSize, random_state=42)
  print('### Dataset split done. Training samples: {}, testing samples: {}'.format(len(toy_train), len(toy_test)))


####### Train & evaluate the models
if trainModels:
  # btagging_models.Fit_RandomForest(toy_train)
  btagging_models.FitKerasModel(toy_train, toy_test, dataname)

minpT = float(sys.argv[3])
maxpT = float(sys.argv[4])

if evalModel:
  btagging_models.EvaluateModelPerformance('Keras_Default', toy_test, dataname, minpT, maxpT)
  # btagging_models.EvaluateModelPerformance('RandomForest', toy_test)


####### Export scores
#if predictData:
#  models = ['Keras_Default', 'RandomForest']
#  for i in range(1,9+1):
#    btagging_helpers.AddScoresToFiles(models, i)
#
#  for i in range(1,20+1):
#    btagging_helpers.AddScoresToFiles_MC(models, i)

if predictData:
  models = ['Keras_Default']
  # btagging_helpers.AddScoresToFiles_MC(models, 1)
  btagging_helpers.AddScoresToFiles(models, 1)
