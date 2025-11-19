#!/usr/bin/env python
from sklearn.model_selection import train_test_split
import btagging_helpers, btagging_models
import pandas as pd
import sys, os
# os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
import tensorflow as tf

#### Settings
recreateData            = False
checkData               = False
trainModels             = True
evalModel               = True
predictData             = False
nClasses                = 3
testSize                = 0.15
valSize                 = 0.15
trainSize               = 0.7
btagging_helpers.gTarget = 2 # 2 = bjets, 1 =cjets

# tf.debugging.set_log_device_placement(True)
# tf.config.set_visible_devices([], 'GPU')

data_MC = pd.DataFrame()
dataname = sys.argv[1]
nJetsSample = int(sys.argv[2])
outputname = sys.argv[5]

if outputname == '':
  outputname = dataname

try:
    os.makedirs("Data/{}".format(dataname), exist_ok=True)
    os.makedirs("Results/{}".format(outputname), exist_ok=True)
    os.makedirs("Models/{}".format(outputname), exist_ok=True)
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
  # toy_train, toy_test = train_test_split(data_MC, test_size=testSize, train_size =trainSize, random_state=42)
  toy_train, hold = train_test_split(data_MC, test_size=(testSize+valSize), train_size =trainSize, random_state=42)  # 70% / 30%
  toy_val, toy_test  = train_test_split(hold,    test_size=0.50,    random_state=42)  # split 30% -> 15%/15%

  print('### Dataset split done. Training samples: {}, testing samples: {}'.format(len(toy_train), len(toy_test)))


####### Train & evaluate the models
if trainModels:
  # btagging_models.Fit_RandomForest(toy_train)
  if nClasses == 2:
    btagging_models.FitKerasModel(toy_train, toy_val, outputname)
  else:
    btagging_models.FitKerasModelMutliClass(toy_train, toy_val, outputname)

minpT = float(sys.argv[3])
maxpT = float(sys.argv[4])

if evalModel:
  if nClasses == 2:
    btagging_models.EvaluateModelPerformance('Keras_Default', toy_test, outputname, minpT, maxpT)
  else:
    btagging_models.EvaluateModelPerformanceMultiClass('Keras_Default', toy_test, outputname, minpT, maxpT)
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
  # btagging_helpers.AddScoresToFiles_MC(models, dataname, outputname)
  btagging_helpers.AddScoresToFiles(models, dataname, outputname)
