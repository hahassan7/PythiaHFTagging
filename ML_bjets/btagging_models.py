from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import RidgeClassifier
import matplotlib.pyplot as plot
import pandas as pd
import btagging_keras, btagging_helpers
import numpy
import tensorflow as tf
from tensorflow.keras.layers import Input
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix

plot.rc('axes', linewidth=2)

#########################################################
def FitModel(label, data, clas, dataname="", sample_weight=None):
  y_train = (data['jetFlavor'] == btagging_helpers.gTarget).values

  X_train = pd.DataFrame(data)
  X_train.drop('jetFlavor', axis=1, inplace=True)

  print('Input parameters: {}'.format(X_train.columns))

  if sample_weight is None:
    clas.fit(X_train, y_train)
  else:
    clas.fit(X_train, y_train, sample_weight = sample_weight)

  # Save the model
  import joblib
  joblib.dump(clas, './Models/{}/Model_{}.pkl'.format(dataname, label))

  print("... done")
  print("")


#########################################################
def Fit_MLP(data):
  print("")
  print("### Training neural network classifier...")
  clas = MLPClassifier(hidden_layer_sizes=(100,100,50), learning_rate='adaptive', verbose=False)
  FitModel('NeuralNetwork_Target{}'.format(btagging_helpers.gTarget), data, clas)

def Fit_RandomForest(data):
  print("")
  print("### Training random forest classifier...")
  clas = RandomForestClassifier(n_estimators=150, max_depth=8)
  FitModel('RandomForest_Target{}'.format(btagging_helpers.gTarget), data, clas)

def Fit_Ridge(data):
  print("")
  print("### Training ridge classifier...")
  clas = RidgeClassifier()
  FitModel('Ridge_Target{}'.format(btagging_helpers.gTarget), data, clas)

def Fit_Keras(data):
  print("")
  print("### Training Keras classifier...")
  FitKerasModel(data, data)


###############################################
def GenerateROCCurve(label, truth, score, suffix = ''):
  """ROC curve & AUC are generated"""
  from sklearn.metrics import roc_curve, roc_auc_score

  currentROCx, currentROCy, _ = roc_curve(truth, score)
  currentAUC = roc_auc_score(truth, score)
  print('AUC={:f}'.format(currentAUC))

  plot.plot(currentROCy, currentROCy, label='Guess ROC', linewidth=2.0)
  plot.plot(currentROCx, currentROCy, label='b-jets', linewidth=2.0)
  plot.subplots_adjust(left=0.15, right=0.95, top=0.90, bottom=0.15)
  plot.legend(loc='lower right', prop={'size': 15})
  plot.tick_params(labelsize=15, length=10, width=2)
  plot.xlim(0, 1)
  plot.ylim(0, 1)
  plot.title('ROC curves', fontsize=20)
  plot.xlabel('False positive rate', fontsize=15)
  plot.ylabel('True positive rate', fontsize=15)
  plot.savefig('./Results/{}/{:s}-ROC-Eval.png'.format(suffix, label), dpi=480)
  plot.clf()


#########################################################
def GetScoreThreshold(scores, efficiency):
  import numpy
  if len(scores.shape) == 2:
    scores = numpy.sort(scores[:, 0])[::-1]
  else:
    scores = numpy.sort(scores)[::-1]
  for i, score in enumerate(scores):
    eff = float(i)/len(scores)
    if eff >= efficiency:
      return score

#########################################################
def GetThresholdEfficiency(scores, threshold):
  import numpy
  if len(scores.shape) == 2:
    scores = numpy.sort(scores[:, 0])[::-1]
  else:
    scores = numpy.sort(scores)[::-1]
  selected = 0
  for score in scores:
    if score >= threshold:
      selected += 1
  return float(selected)/len(scores)


#########################################################
def EvaluateModelPerformance(model, data_in, name, minpT, maxpT):
  print('\n###### Evaluating {}'.format(model))

  ### Filter + split datasets
  data = pd.DataFrame(data_in.loc[(data_in['jetPt'] >= minpT) & (data_in['jetPt'] < maxpT)])
  data_test_B  = pd.DataFrame(data.loc[(data['jetFlavor'] == 2)])
  data_test_C  = pd.DataFrame(data.loc[(data['jetFlavor'] == 1)])
  data_test_LF = pd.DataFrame(data.loc[(data['jetFlavor'] == 0)])

  ### Combined dataset (with equal class weights)
  data_combined = pd.concat([data_test_B, data_test_C, data_test_LF])
  truth_test_combined = (data_combined['jetFlavor'] == btagging_helpers.gTarget).values
  data_test_combined = pd.DataFrame(data_combined)

  data_test_B.drop('jetFlavor', axis=1, inplace=True)
  data_test_C.drop('jetFlavor', axis=1, inplace=True)
  data_test_LF.drop('jetFlavor', axis=1, inplace=True)
  data_test_combined.drop('jetFlavor', axis=1, inplace=True)

  ### Extract scores for sklearn classifier
  if not 'keras' in model.lower():
    import joblib
    clas = joblib.load('./Models/{}/Model_{}.pkl'.format(name, model))

    score_combined = clas.predict_proba(data_test_combined)[:,1]
    score_B  = clas.predict_proba(data_test_B)[:,1]
    score_C  = clas.predict_proba(data_test_C)[:,1]
    score_LF = clas.predict_proba(data_test_LF)[:,1]
    y_pred_combined = clas.predict(data_test_combined)

  ### Extract scores for keras classifier
  else:
    X_test_combined  = GetDataForKeras(data_test_combined)
    X_test_B  = GetDataForKeras(data_test_B)
    X_test_C  = GetDataForKeras(data_test_C)
    X_test_LF = GetDataForKeras(data_test_LF)
    myModel = btagging_keras.AliMLKerasModel(2)
    myModel.LoadModel('{}/Model_{}'.format(name, model))#_Default')
    score_combined  = myModel.fModel.predict(X_test_combined, batch_size=512, verbose=0)
    score_B  = myModel.fModel.predict(X_test_B, batch_size=512, verbose=0)
    score_C  = myModel.fModel.predict(X_test_C, batch_size=512, verbose=0)
    score_LF = myModel.fModel.predict(X_test_LF, batch_size=512, verbose=0)

  threshold = GetScoreThreshold(score_B, 0.6)
  print('\n##########################')
  print('\nThresold score is: {} \n'.format(threshold))
  print('b efficiency: {}'.format(GetThresholdEfficiency(score_B, threshold)))
  print('c efficiency: {}'.format(GetThresholdEfficiency(score_C, threshold)))
  print('udsg efficiency: {}'.format(GetThresholdEfficiency(score_LF, threshold)))
  print('##########################')

  # Evaluate model
  GenerateROCCurve(model, truth_test_combined, score_combined, name)
  SaveHistogram('./Results/{}/{:s}-Scores.png'.format(name, model), 'Scores on validation data', tuple([score_B, score_C, score_LF]), ('B','C', "light flavor"), rangex=(0,1), logY=True)

  # Save confusion matrix
  y_pred_combined = (score_combined >= threshold).astype(int)
  SaveConfusionMatrix(truth_test_combined, y_pred_combined, classes=['Light/Charm jets', 'Beauty jets'], fname='./Results/{}/{:s}-ConfusionMatrix.png'.format(name, model))


#########################################################
def SaveHistogram(fname, title, y, functionlabels, rangex=(0,1), legendloc='upper right', logY=False, nbins=100):
  import math
  """Simple histogram helper function"""
  # Check input data
  for i in range(len(y)):
    for item in y[i]:
      if math.isnan(item):
        print('Histogram {:s} was not saved due to invalid values.'.format(fname))
        return

  for i in range(len(y)):
    plot.hist(y[i], label=functionlabels[i], bins=nbins, histtype='step', density=True, linewidth=2.0)
  plot.subplots_adjust(left=0.15, right=0.95, top=0.90, bottom=0.15)
  plot.tick_params(labelsize=15, length=10, width=2)
  plot.legend(loc=legendloc, prop={'size': 15})
  plot.title(title, fontsize=20, pad=10)
  plot.xlabel('Scores', fontsize=15)
  plot.xlim(rangex[0], rangex[1])
#  plot.yscale('log', nonposy='clip')

  if logY:
    plot.yscale('log')

  plot.savefig("{:s}".format(fname), dpi=480)
  #plot.show()
  plot.clf()

#########################################################
def SaveConfusionMatrix(y_true, y_pred, classes, fname, normalize=False, title='Confusion matrix', cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = confusion_matrix(y_true, y_pred, normalize='true')
    else:
        cm = confusion_matrix(y_true, y_pred)

    plt.figure(figsize=(10, 7))
    sns.heatmap(cm, annot=True, fmt='.2f' if normalize else 'd', cmap=cmap, xticklabels=classes, yticklabels=classes)
    plt.title(title)
    plt.xlabel('Predicted label')
    plt.ylabel('True label')
    plt.savefig(fname, dpi=480)
    plt.close()

#########################################################
def FitKerasModel(data, val_data, dataname="", verbose=True):

  ### Prepare data
  X_train = GetDataForKeras(data)
  X_test = GetDataForKeras(val_data)
  y_train = (data['jetFlavor'] == btagging_helpers.gTarget).values
  print("tain shape is: {}".format(y_train.shape))
  print(y_train)
  y_test = (val_data['jetFlavor'] == btagging_helpers.gTarget).values

  ### Prepare model
  myModel = btagging_keras.AliMLKerasModel(2)
  inputLayer = Input(shape=(6,), name='B{:d}_GV_In'.format(len(myModel.fTempModelBranches)))
  myModel.fModelBranchesOutput.append('B{:d}_GV_{:d}_{:3.2f}_Out'.format(len(myModel.fTempModelBranches), 1, 0))
  myModel.fModelBranchesInput.append(inputLayer.name)
  myModel.fTempModelBranchesInput.append(inputLayer)
  myModel.fTempModelBranches.append(inputLayer)

  myModel.AddBranchCNN1D([64,32,32,16], [2,0,0,0], 1, [11,7,2,2], 0.2, inputShape=(10, 9))
  myModel.AddBranchCNN1D([64,32,32,16], [2,0,0,0], 1, [11,7,2,2], 0.2, inputShape=(10, 10))
  myModel.SetFinalLayerVariable([200,200,200,100],0.2, 'ridge')

  myModel.fInit = 'he_uniform'
  myModel.fOptimizer = 'adam'
  myModel.fLossFunction = 'binary_crossentropy'
  myModel.fBatchSize = 1000
  myModel.CreateModel('{}/Model_Keras_Default'.format(dataname))
  myModel.fLearningRate = 0.001
  myModel.PrintProperties()

  ### Train model
  myModel.TrainModel(X_train, y_train, X_test, y_test, numEpochs = 10)
  myModel.SaveModel()


#########################################################
def GetDataForKeras(data):


  arr_branch0 = numpy.array([data['jetPt'], data['jetEta'], data['jetPhi'], data['nTracks'], data['nSV'], data['jetMass']])
  print("Array0 shape: {}".format(arr_branch0.shape) )
  arr_branch0 = numpy.swapaxes(arr_branch0, 0, 1)
  print("Array0 shape after: {}".format(arr_branch0.shape) )

  arr_mTrackpT        = pd.concat([data['trackPt[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mTrackEta   = pd.concat([data['trackEta[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mDotProdTrackJet = pd.concat([data['dotProdTrackJet[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mDotProdTrackJetOverJet = pd.concat([data['dotProdTrackJetOverJet[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mDeltaRJetTrack = pd.concat([data['deltaRJetTrack[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mSignedIP2D = pd.concat([data['signedIP2D[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mSignedIP3D = pd.concat([data['signedIP3D[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mMomFraction = pd.concat([data['momFraction[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mDeltaRTrackVertex = pd.concat([data['deltaRTrackVertex[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_branch1 = numpy.array([arr_mTrackpT, arr_mTrackEta,arr_mDotProdTrackJet, arr_mDotProdTrackJetOverJet, arr_mDeltaRJetTrack, arr_mSignedIP2D, arr_mSignedIP3D, arr_mMomFraction, arr_mDeltaRTrackVertex])
  print("Array1 shape: {}".format(arr_branch1.shape) )
  arr_branch1 = numpy.swapaxes(arr_branch1, 0, 1)
  print("Array1 shape: {}".format(arr_branch1.shape) )
  arr_branch1 = numpy.swapaxes(arr_branch1, 1, 2)
  print("Array1 shape: {}".format(arr_branch1.shape) )

  arr_mSVpT        = pd.concat([data['svPt[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mDeltaRSVJet     = pd.concat([data['deltaRSVJet[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mSVMass     = pd.concat([data['svMass[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mSVfE     = pd.concat([data['svfE[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mIPXY     = pd.concat([data['svIPxy[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mCPA     = pd.concat([data['svCPA[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mChi2PCA     = pd.concat([data['svChi2PCA[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mDispersion     = pd.concat([data['svDispersion[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mDecayLength2D     = pd.concat([data['svDecayLength2D[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_mDecayLength3D     = pd.concat([data['svDecayLength3D[{}]'.format(i)] for i in range(10)], axis=1).values
  arr_branch2 = numpy.array([arr_mSVpT, arr_mDeltaRSVJet, arr_mSVMass, arr_mSVfE, arr_mIPXY, arr_mCPA, arr_mChi2PCA, arr_mDispersion, arr_mDecayLength2D, arr_mDecayLength3D])
  print("Array2 shape: {}".format(arr_branch2.shape) )
  arr_branch2 = numpy.swapaxes(arr_branch2, 0, 1)
  print("Array2 shape: {}".format(arr_branch2.shape) )
  arr_branch2 = numpy.swapaxes(arr_branch2, 1, 2)
  print("Array2 shape: {}".format(arr_branch2.shape) )

  return [arr_branch0,arr_branch1,arr_branch2]
