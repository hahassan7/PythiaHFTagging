import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers import concatenate, LSTM, Dense, Dropout, Activation, Conv1D, MaxPooling1D, Flatten, Conv2D, MaxPooling2D, SimpleRNN, ZeroPadding2D, Input
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.callbacks import History, EarlyStopping
from tensorflow.keras.regularizers import l1_l2, l1, l2
from tensorflow.keras.utils import plot_model
import tensorflow.keras.backend as K
import matplotlib.pyplot as plot

plot.rc('axes', linewidth=3)

import sys, os, math, numpy, pickle, logging, ROOT

#__________________________________________________________________________________________________________
class AliMLKerasModel:
  """Keras model meta class"""
  ###############################################
  def __init__(self, numClasses):
    self.fModel = None
    # Temp branches used in building the model
    self.fTempModelBranches = []
    self.fTempModelBranchesInput  = []
    # Humanreadable properties of the network
    self.fModelBranchesInput  = []
    self.fModelBranchesOutput = []
    # Final layer properties (top vanilla network)
    self.fFinalLayerStructure = []
    self.fFinalLayerNumLayers = 3
    self.fFinalLayerNumNeuronsPerLayer = 512
    self.fFinalLayerDropout = 0.5
    self.fFinalLayerRegularization = 'elasticnet'
    self.fFinalLayerActivation = 'relu'
    # Model properties
    self.fOptimizer = 'SGD'
    self.fLossFunction = 'binary_crossentropy'
    self.fModelName = 'NonameMetaModel'
    self.fLearningRate = 0.05
    self.fInit = 'he_normal'
    self.fNumClasses = numClasses
    self.fBatchSize = 512
    # Misc
    self.fPerEpochValidation = True # do validation per epoch
    self.fShowModelSummary = True
    self.fRequestedData = []
    self.fResults = None

  ###############################################
  def PrintProperties(self):
    epochs = len(self.fResults.fHistoryLearningRate)
    # General information on the model
    if epochs > 0:
      if self.fNumClasses == 2: # bin classification task
          print('\nLoss={:4.3f}, Acc={:4.3f}, AUC={:4.3f} after {:3d} epochs (lr={:6.5f}). optimizer={:s}, loss function={:s}, init={:s}'.format(self.fResults.GetReachedLoss(), self.fResults.GetReachedAccuracy(), self.fResults.GetReachedAUC(), epochs, self.fResults.fHistoryLearningRate[epochs-1], self.fOptimizer, self.fLossFunction, self.fInit))
      else: # classification task/regression task
          print('\nLoss={:4.3f} after {:3d} epochs (lr={:6.5f}). optimizer={:s}, loss function={:s}, init={:s}'.format(self.fResults.GetReachedLoss(), epochs, self.fResults.fHistoryLearningRate[epochs-1], self.fOptimizer, self.fLossFunction, self.fInit))

    print('\n###########################\nUsing {:s}. Model branches:'.format(self.fModelName))
    for inBranch in self.fModelBranchesInput:
      print('  {}'.format(inBranch))
    print('###########################\n')

  ###############################################
  def GetOutputLayer(self, numClasses):
    if numClasses == 1: # regression task
      return Dense(1, activation='linear')
    elif numClasses == 2: # bin classification task
      return Dense(1, activation='sigmoid', kernel_initializer='he_normal')
    elif numClasses > 2: # classification task
      return Dense(numClasses, activation='sigmoid', kernel_initializer='he_normal')

  ###############################################
  def GetResultsObject(self, numClasses, name):
    if numClasses == 1:
      results = AliMLModelResultsRegressionTask(name)
    elif numClasses == 2:
      results = AliMLModelResultsBinClassifier(name)
    elif numClasses > 2:
      results = AliMLModelResultsMultiClassifier(name, numClasses)

    return results

  ###############################################
  def CreateModel(self, mname):
    if self.fModel:
      raise ValueError('Model already exists. Create a new one instead.')

    if not len(self.fTempModelBranches):
      raise ValueError('No branches given. Add some branches.')

    self.fModelName = mname

    # Merge the model branches to one model (that's the output we want to train)
    if len(self.fTempModelBranches) > 1:
      modelOutput = concatenate([self.fTempModelBranches[i] for i in range(len(self.fTempModelBranches))], name='ConcatenateLayer')
    else:
      modelOutput = self.fTempModelBranches[0]

    if self.fFinalLayerRegularization == 'elasticnet':
      regul = l1_l2(0.001)
    elif self.fFinalLayerRegularization == 'ridge':
      regul = l2(0.001)
    elif self.fFinalLayerRegularization == 'lasso':
      regul = l1(0.001)
    else:
      raise ValueError('Regularization mode {:s} not recognized'.format(self.fFinalLayerRegularization))

    # Add final vanilla dense layer
    if self.fFinalLayerStructure == []:
      for i in range(self.fFinalLayerNumLayers):
        modelOutput = Dense(self.fFinalLayerNumNeuronsPerLayer, activation=self.fFinalLayerActivation, kernel_initializer=self.fInit, kernel_regularizer=regul, name='Final_{:d}_{:d}_activation_{:s}_regularization_{:s}'.format(i, self.fFinalLayerNumNeuronsPerLayer, self.fFinalLayerActivation, self.fFinalLayerRegularization))(modelOutput)

        if self.fFinalLayerDropout:
          modelOutput = Dropout(self.fFinalLayerDropout, name='Final_{:d}_{:3.2f}'.format(i, self.fFinalLayerDropout))(modelOutput)
    else:
      for i, nNodes in enumerate(self.fFinalLayerStructure):
        modelOutput = Dense(nNodes, activation=self.fFinalLayerActivation, kernel_initializer=self.fInit, kernel_regularizer=regul, name='Final_{:d}_{:d}_activation_{:s}_regularization_{:s}'.format(i, self.fFinalLayerNumNeuronsPerLayer, self.fFinalLayerActivation, self.fFinalLayerRegularization))(modelOutput)

        if self.fFinalLayerDropout:
          modelOutput = Dropout(self.fFinalLayerDropout, name='Final_{:d}_{:3.2f}'.format(i, self.fFinalLayerDropout))(modelOutput)


    # Add final output layer
    modelOutput = self.GetOutputLayer(self.fNumClasses)(modelOutput)
    # Create model proxy
    self.fModel = Model(inputs=[inBranch for inBranch in self.fTempModelBranchesInput], outputs=modelOutput)

    ### *COMPILE*
    self.fModel.compile(loss=self.fLossFunction, optimizer=self.fOptimizer, metrics=["accuracy"])

    if self.fShowModelSummary:
      self.fModel.summary()
      plot_model(self.fModel, to_file='./Models/{:s}.png'.format(mname))

    self.fResults = self.GetResultsObject(self.fNumClasses, self.fModelName)

    return self.fModel

  ###############################################
  def TrainModel(self, data, truth, validationData, validationTruth, model=None, results=None, numEpochs=1):
    callbacks = []
    callbacks.append(AliMLKerasModel_EnhancedProgbarLogger())
    # Early stopping
    #callbacks.append(EarlyStopping(monitor='val_loss', patience=4))
    # Learning rate reduction on plateau
    callbacks.append(keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.1, patience=20, verbose=0))

    if model == None:
      model = self.fModel

    if results == None:
      results = self.fResults

    callbacks.append(AliMLKerasModel_CallbackSaveModel(self.SaveModel, model, validationData, validationTruth, self.fBatchSize, results, len(data[0]), self.fPerEpochValidation))

    ###############
    # Set learning should work for theano & tensorflow
    K.set_value(model.optimizer.learning_rate.value, self.fLearningRate)

    # Train and test for numEpochs epochs
    for i in range(len(data)):
      print("shape of data: {}".format(data[i].shape))
    hist = (model.fit(data, truth, epochs=numEpochs, batch_size=self.fBatchSize,
                  validation_data=(validationData, validationTruth), verbose=0, callbacks=callbacks))

  ###############################################
  def SaveModel(self):
    saveObj = {}
    # Prevent saving of the keras model with pickle
    # Save keras model separately
    for key in self.__dict__:
      if key == 'fModel':
        if self.fModel:
          self.fModel.save('./Models/{:s}.h5'.format(self.fModelName))
      elif key.startswith('fTemp'):
        pass
      elif key == 'fOptimizer' and not isinstance(self.fOptimizer, str):
        saveObj[key] = self.fOptimizer.__class__.__name__
      else:
        saveObj[key] = self.__dict__[key]

    if sys.version_info <= (3, 12):
      saved_model_path = "./Models/{:s}/saved_model".format(self.fModelName.split('/')[0])
      tf.saved_model.save(self.fModel, saved_model_path)

    pickle.dump(saveObj, open('./Models/{:s}.p'.format(self.fModelName), 'wb'))

  ###############################################
  def LoadModel(self, fname, meta_data_only=False):
    if not os.path.isfile('./Models/{:s}.p'.format(fname)) or not os.path.isfile('./Models/{:s}.h5'.format(fname)):
       raise ValueError('Error: Cannot load model {:s} due to missing files!'.format(fname))

    self.__dict__ = pickle.load(open('./Models/{:s}.p'.format(fname), 'rb'))
    if not meta_data_only:
      self.fModel = load_model('./Models/{:s}.h5'.format(fname))
    self.fModelName = fname

  ###############################################
  def ExportModel(self):
    """Export model to json & weights-only h5 file. Usable in LWTNN"""
    self.fModel.save_weights('./Models/{:s}.weights.h5'.format(self.fModelName), overwrite=True)

    # Save network architecture to file
    jsonfile = open('./Models/{:s}.json', 'w')
    jsonfile.write(self.fModel.to_json())
    jsonfile.close()


  ###############################################
  def AddBranchDense(self, nLayers, nNeuronsPerLayer, dropout, inputShape, activation='relu', regularization='off'):
    self.fRequestedData.append(inputShape)
    branchID = len(self.fTempModelBranches)

    if regularization   == 'elasticnet':
      regul = l1_l2(0.001)
    elif regularization == 'ridge':
      regul = l2(0.001)
    elif regularization == 'lasso':
      regul = l1(0.001)
    elif regularization == 'off':
      regul = None
    else:
      raise ValueError('Regularization mode {:s} not recognized'.format(regul))

    # Create fully-connected layers
    inputLayer = Input(shape=inputShape, name='B{:d}_FC'.format(branchID))
    for i in range(nLayers):
      if i==0:
        model = Dense(nNeuronsPerLayer, activation=activation, kernel_initializer=self.fInit, kernel_regularizer=regul, name='B{:d}_FC_{:d}_{:d}_activation_{:s}_shape_{:d}'.format(branchID, i, nNeuronsPerLayer, activation, inputShape))(inputLayer)
      else:
        model = Dense(nNeuronsPerLayer, activation=activation, kernel_initializer=self.fInit, kernel_regularizer=regul, name='B{:d}_FC_{:d}_{:d}_activation_{:s}'.format(branchID, i, nNeuronsPerLayer, activation))(model)

      if dropout:
        model = Dropout(dropout, name='B{:d}_FC_{:d}_{:3.2f}'.format(branchID, i, dropout))(model)

    self.fModelBranchesOutput.append('B{:d}_FC_{:d}_{:3.2f}'.format(branchID, i, dropout))
    self.fModelBranchesInput.append(inputLayer.name)
    self.fTempModelBranchesInput.append(inputLayer)
    self.fTempModelBranches.append(model)

  ###############################################
  def AddBranchDenseCustom(self, layers, dropout, inputShape, activation='relu', regularization='off'):
    self.fRequestedData.append(inputShape)
    branchID = len(self.TempfModelBranches)

    if regularization   == 'elasticnet':
      regul = l1_l2(0.001)
    elif regularization == 'ridge':
      regul = l2(0.001)
    elif regularization == 'lasso':
      regul = l1(0.001)
    elif regularization == 'off':
      regul = None
    else:
      raise ValueError('Regularization mode {:s} not recognized'.format(regul))

    # Create fully-connected layers
    inputLayer = Input(shape=inputShape, name='B{:d}_FC'.format(branchID))
    for i in range(len(layers)):
      if i==0:
        model = Dense(layers[i], activation=activation, kernel_initializer=self.fInit, kernel_regularizer=regul, name='B{:d}_FC_{:d}_{:d}_activation_{:s}'.format(branchID, i, layers[i], activation))(inputLayer)
      else:
        model = Dense(layers[i], activation=activation, kernel_initializer=self.fInit, kernel_regularizer=regul, name='B{:d}_FC_{:d}_{:d}_activation_{:s}'.format(branchID, i, layers[i], activation))(model)

      if dropout:
        model = Dropout(dropout, name='B{:d}_FC_{:d}_{:3.2f}'.format(branchID, i, dropout))(model)

    self.fModelBranchesOutput.append('B{:d}_FC_{:d}_{:3.2f}'.format(branchID, i, dropout))
    self.fModelBranchesInput.append(inputLayer.name)
    self.fTempModelBranchesInput.append(inputLayer)
    self.fTempModelBranches.append(model)

  ###############################################
  def AddBranchCNN1D(self, seqConvFilters, seqMaxPoolings, subsampling, seqKernelSizes, dropout, inputShape, activation='relu'):
    self.fRequestedData.append(inputShape)
    branchID = len(self.fTempModelBranches)

    inputLayer = Input(shape=inputShape, name='B{:d}_CNN1D'.format(branchID))
    for i in range(len(seqConvFilters)):
      if i==0:
        model = Conv1D(seqConvFilters[i], seqKernelSizes[i], activation=activation, kernel_initializer=self.fInit, strides=subsampling, padding='same', name=self.GetCompatibleName('B{:d}_CNN1D_{:d}_Kernels{}_Stride_{:d}_activation_{:s}'.format(branchID, i, seqKernelSizes, subsampling, activation)))(inputLayer)
      else:
        model = Conv1D(seqConvFilters[i], seqKernelSizes[i], activation=activation, kernel_initializer=self.fInit, padding='same', name=self.GetCompatibleName('B{:d}_CNN1D_{:d}_Kernels{}_activation_{:s}'.format(branchID, i, seqKernelSizes, activation)))(model)

      if seqMaxPoolings[i] > 0:
        model = MaxPooling1D(pool_size=seqMaxPoolings[i], name='B{:d}_CNN1D_{:d}_{}'.format(branchID, i, seqMaxPoolings[i]))(model)
      if dropout:
        model = Dropout(dropout, name='B{:d}_CNN1D_{:d}_{:3.2f}'.format(branchID, i, dropout))(model)
    model = Flatten(name='B{:d}_CNN1D_Output'.format(branchID))(model)

    self.fModelBranchesOutput.append('B{:d}_CNN1D_Output'.format(branchID))
    self.fModelBranchesInput.append(inputLayer.name)
    self.fTempModelBranchesInput.append(inputLayer)
    self.fTempModelBranches.append(model)

  ###############################################
  def AddBranchLC1D(self, seqConvFilters, seqMaxPoolings, subsampling, seqKernelSizes, dropout, inputShape, activation='relu'):
    self.fRequestedData.append(inputShape)
    branchID = len(self.fTempModelBranches)

    inputLayer = Input(shape=inputShape, name='B{:d}_LC1D'.format(branchID))
    for i in range(0, len(seqConvFilters)):
      if i==0:
        model = Dense(seqConvFilters[i], seqKernelSizes[i], activation=activation, kernel_initializer=self.fInit, strides=subsampling, padding='valid', name=self.GetCompatibleName('B{:d}_LC1D_{:d}_Kernels{}_Stride_{:d}_activation_{:s}'.format(branchID, i, seqKernelSizes, subsampling, activation)))(inputLayer)
      else:
        model = Dense(seqConvFilters[i], seqKernelSizes[i], activation=activation, kernel_initializer=self.fInit, padding='valid', name=self.GetCompatibleName('B{:d}_LC1D_{:d}_Kernels{}_activation_{:s}'.format(branchID, i, seqKernelSizes, activation)))(model)

      if seqMaxPoolings[i] > 0:
        model = MaxPooling1D(pool_size=seqMaxPoolings[i], name='B{:d}_LC1D_{:d}_{}'.format(branchID, i, seqMaxPoolings[i]))(model)
      if dropout:
        model = Dropout(dropout, name='B{:d}_LC1D_{:d}_{:3.2f}'.format(branchID, i, dropout))(model)
    model = Flatten(name='B{:d}_LC1D_Output'.format(branchID))(model)

    self.fModelBranchesOutput.append('B{:d}_LC1D_Output'.format(branchID))
    self.fModelBranchesInput.append(inputLayer.name)
    self.fTempModelBranchesInput.append(inputLayer)
    self.fTempModelBranches.append(model)

  ###############################################
  def AddBranchCNN2D(self, seqConvFilters, seqMaxPoolings, subsampling, seqKernelSizes, dropout, inputShape, activation='relu'):
    self.fRequestedData.append(inputShape)
    branchID = len(self.fTempModelBranches)

    inputLayer = Input(shape=inputShape, name='B{:d}_CNN2D'.format(branchID))
    for i in range(0, len(seqConvFilters)):
      if i==0:
        model = Conv2D(seqConvFilters[i], seqKernelSizes[i], activation=activation, kernel_initializer=self.fInit, strides=subsampling, padding='same', name=self.GetCompatibleName('B{:d}_CNN2D_{:d}_Kernels{}_Stride_{:d}_activation_{:s}'.format(branchID, i, seqKernelSizes, subsampling, activation)))(inputLayer)
      else:
        model = Conv2D(seqConvFilters[i], seqKernelSizes[i], activation=activation, kernel_initializer=self.fInit, padding='same', name=self.GetCompatibleName('B{:d}_CNN2D_{:d}_Kernels{}_activation_{:s}'.format(branchID, i, seqKernelSizes, activation)))(model)

      if seqMaxPoolings[i] > 0:
        model = MaxPooling2D(pool_size=(seqMaxPoolings[i], seqMaxPoolings[i]), name='B{:d}_CNN2D_{:d}_{}'.format(branchID, i, seqMaxPoolings[i]))(model)
      if dropout:
        model = Dropout(dropout, name='B{:d}_CNN2D_{:d}_{:3.2f}'.format(branchID, i, dropout))(model)
    model = Flatten(name='B{:d}_CNN2D_Output'.format(branchID))(model)

    self.fModelBranchesOutput.append('B{:d}_CNN2D_Output'.format(branchID))
    self.fModelBranchesInput.append(inputLayer.name)
    self.fTempModelBranchesInput.append(inputLayer)
    self.fTempModelBranches.append(model)

  ###############################################
  def AddBranchLC2D(self, seqConvFilters, seqMaxPoolings, subsampling, seqKernelSizes, dropout, inputShape, activation='relu'):
    self.fRequestedData.append(inputShape)
    branchID = len(self.fTempModelBranches)

    inputLayer = Input(shape=inputShape, name='B{:d}_LC2D'.format(branchID))
    for i in range(0, len(seqConvFilters)):
      if i==0:
        model = LocallyConnected2D(seqConvFilters[i], seqKernelSizes[i], activation=activation, kernel_initializer=self.fInit, strides=subsampling, padding='valid', name=self.GetCompatibleName('B{:d}_LC2D_{:d}_Kernels{}_Stride_{:d}_activation_{:s}'.format(branchID, i, seqKernelSizes, subsampling, activation)))(inputLayer)
      else:
        model = LocallyConnected2D(seqConvFilters[i], seqKernelSizes[i], activation=activation, kernel_initializer=self.fInit, padding='valid', name=self.GetCompatibleName('B{:d}_LC2D_{:d}_Kernels{}_activation_{:s}'.format(branchID, i, seqKernelSizes, activation)))(model)

      if seqMaxPoolings[i] > 0:
        model = MaxPooling2D(pool_size=(seqMaxPoolings[i], seqMaxPoolings[i]), name='B{:d}_LC2D_{:d}_{}'.format(branchID, i, seqMaxPoolings[i]))(model)
      if dropout:
        model = Dropout(dropout, name='B{:d}_LC2D_{:d}_{:3.2f}'.format(branchID, i, dropout))(model)
    model = Flatten(name='B{:d}_LC2D_Output'.format(branchID))(model)

    self.fModelBranchesOutput.append('B{:d}_LC2D_Output'.format(branchID))
    self.fModelBranchesInput.append(inputLayer.name)
    self.fTempModelBranchesInput.append(inputLayer)
    self.fTempModelBranches.append(model)


  ###############################################
  def AddBranchLSTM(self, outputDim, nLayers, dropout, inputShape, activation='tanh'):
    self.fRequestedData.append(inputShape)
    branchID = len(self.fTempModelBranches)

    inputLayer = Input(shape=inputShape, name='B{:d}_LSTM'.format(branchID))
    for i in range(nLayers):
      if i==nLayers-1: # last layer
        model = LSTM(outputDim, activation=activation, name='B{:d}_LSTM_{:d}_{:d}_activation_{:s}'.format(branchID, i, outputDim, activation))(inputLayer if i==0 else model)
      else: # other layers
        model = LSTM(outputDim, return_sequences=True, activation=activation, name='B{:d}_LSTM_{:d}_{:d}_activation_returnSeq_{:s}'.format(branchID, i, outputDim, activation))(inputLayer if i==0 else model)
      
      if dropout:
        model = Dropout(dropout, name='B{:d}_LSTM_{:d}_{:3.2f}'.format(branchID, i, dropout))(model)

    self.fModelBranchesOutput.append('B{:d}_LSTM_{:d}_{:3.2f}'.format(branchID, i, dropout))
    self.fModelBranchesInput.append(inputLayer.name)
    self.fTempModelBranchesInput.append(inputLayer)
    self.fTempModelBranches.append(model)

  ###############################################
  def AddBranchRNN(self, outputDim, nLayers, dropout, inputShape, activation='tanh'):
    self.fRequestedData.append(inputShape)
    branchID = len(self.fTempModelBranches)

    inputLayer = Input(shape=inputShape, name='B{:d}_LSTM'.format(branchID))
    for i in range(nLayers):
      if i==nLayers-1: # last layer
        model = SimpleRNN(outputDim, activation=activation, name='B{:d}_LSTM_{:d}_{:d}_activation_{:s}'.format(branchID, i, outputDim, activation))(inputLayer if i==0 else model)
      else: # other layers
        model = SimpleRNN(outputDim, return_sequences=True, activation=activation, name='B{:d}_LSTM_{:d}_{:d}_activation_returnSeq_{:s}'.format(branchID, i, outputDim, activation))(inputLayer if i==0 else model)

      if dropout:
        model = Dropout(dropout, name='B{:d}_LSTM_{:d}_{:3.2f}'.format(branchID, i, dropout))(model)

    self.fModelBranchesOutput.append('B{:d}_LSTM_{:d}_{:3.2f}'.format(branchID, i, dropout))
    self.fModelBranchesInput.append(inputLayer.name)
    self.fTempModelBranchesInput.append(inputLayer)
    self.fTempModelBranches.append(model)

  ###############################################
  def SetFinalLayer(self, nLayers, nNeuronsPerLayer, dropout, regularization='ridge', activation='relu'):
    self.fFinalLayerStructure = []
    self.fFinalLayerNumLayers = nLayers
    self.fFinalLayerNumNeuronsPerLayer = nNeuronsPerLayer
    self.fFinalLayerDropout = dropout
    self.fFinalLayerRegularization = regularization
    self.fFinalLayerActivation = activation

  ###############################################
  def SetFinalLayerVariable(self, layers, dropout, regularization='ridge', activation='relu'):
    self.fFinalLayerStructure = layers
    self.fFinalLayerNumLayers = 0
    self.fFinalLayerNumNeuronsPerLayer = 0
    self.fFinalLayerDropout = dropout
    self.fFinalLayerRegularization = regularization
    self.fFinalLayerActivation = activation

  ###############################################
  def GetCompatibleName(self, name):
    return name.replace(', ', '_').replace('[','_').replace(']_','_').replace(']','_').replace(' ','_')

#__________________________________________________________________________________________________________
#__________________________________________________________________________________________________________
class AliMLKerasModel_CallbackSaveModel(keras.callbacks.Callback):
  """Keras callback that call the SaveModel function after each epoch"""
  def __init__(self, saveFunction, model, valData, valTruth, batch_size, results, nevents, dovalidation):
    super(AliMLKerasModel_CallbackSaveModel, self).__init__()
    self.fSaveFunction = saveFunction
    self.fModel = model
    self.fValidationData  = valData
    self.fValidationTruth = valTruth
    self.fBatchSize = batch_size
    self.fResults = results
    self.fNEvents = nevents
    self.fDoValidation = dovalidation
  def on_epoch_end(self, epoch, logs=None):
    # Save the model on each epoch's end
    self.fSaveFunction()
    if self.fDoValidation:
      learningRate = K.get_value(self.fModel.optimizer.learning_rate.value)
      self.fResults.AddResult(logs.get('loss'), logs.get('val_loss'), learningRate, self.fNEvents, logs.get('acc'), logs.get('val_acc'), self.fModel, self.fValidationData, self.fValidationTruth, self.fBatchSize)

class AliMLKerasModel_EnhancedProgbarLogger(keras.callbacks.ProgbarLogger):
  """Custom progress bar"""
  def on_train_begin(self, logs=None):
    keras.callbacks.ProgbarLogger.on_train_begin(self, logs)
    self.verbose = 1
  def on_epoch_begin(self, epoch, logs=None):
    print('')
    logging.info('')
    keras.callbacks.ProgbarLogger.on_epoch_begin(self, epoch, logs)

#__________________________________________________________________________________________________________




class AliMLModelResultsBase:
  """Base class that holds model results and perform tests on it"""
  ###############################################
  def __init__(self, name):
    self.fModelName = name

    # Arrays holding values obtained for each learning step
    self.fHistoryLossTraining       = []
    self.fHistoryLossValidation     = []
    self.fHistoryLearningRate       = []
    self.fHistoryNumberEvents       = []


  ###############################################
  def GetReachedLoss(self):
    return self.fHistoryLossValidation[len(self.fHistoryLossValidation)-1]


  ###############################################
  def AddResult(self, loss_train, loss_val, lr, numEvents):
    """Add performance results from a model training"""
    self.fHistoryLossTraining.append(loss_train)
    self.fHistoryLossValidation.append(loss_val)
    self.fHistoryLearningRate.append(lr)

    # In case, the result contains more than one epoch, create the list of numEvents
    offset = 0 if len(self.fHistoryNumberEvents) == 0 else self.fHistoryNumberEvents[len(self.fHistoryNumberEvents)-1]
    self.fHistoryNumberEvents.append(offset+numEvents)


  ###############################################
  def CreatePlots(self):
    """This function saves the plots using the existing results data"""
    ##### Save histograms/plots
    epochs    = range(len(self.fHistoryNumberEvents))
    # Loss
    SavePlot('./Results/{:s}-Loss.png'.format(self.fModelName), 'Loss function', x=epochs, y=(self.fHistoryLossTraining, self.fHistoryLossValidation), rangey=(0.8* min(self.fHistoryLossValidation),1.2*max(self.fHistoryLossValidation)), functionlabels=('Training', 'Validation'), axislabels=('Epoch', 'Loss'))


#__________________________________________________________________________________________________________
class AliMLModelResultsBinClassifier(AliMLModelResultsBase):
  """Model result class for binary classifier results"""

  ###############################################
  def __init__(self, name):
    AliMLModelResultsBase.__init__(self, name)

    self.fCurrentAUC                = None
    self.fCurrentROCx               = None
    self.fCurrentROCy               = None
    self.fCurrentTestScores         = None
    # Values for classification
    self.fHistoryAccuracyTraining   = []
    self.fHistoryAccuracyValidation = []
    self.fHistoryAUC                = []


  ###############################################
  def GetReachedAccuracy(self):
    return self.fHistoryAccuracyValidation[len(self.fHistoryAccuracyValidation)-1]
  def GetReachedAUC(self):
    return self.fHistoryAUC[len(self.fHistoryAUC)-1]


  ###############################################
  def AddResult(self, loss_train, loss_val, lr, numEvents, accuracy_train, accuracy_val, model, data, truth, test_batch_size):
    """Add performance results from a model training"""
    AliMLModelResultsBase.AddResult(self, loss_train, loss_val, lr, numEvents)
    self.fHistoryAccuracyTraining.append(accuracy_train)
    self.fHistoryAccuracyValidation.append(accuracy_val)

    if not model:
      raise ValueError('Cannot test a model that was not correctly created.')


    self.RunTests(model, data, truth, test_batch_size)
    self.fHistoryAUC.append(self.fCurrentAUC)

    self.CreatePlots()


  ###############################################
  def RunTests(self, model, data, truth, test_batch_size):
    """Runs tests"""

    # Calculate ROC/ROC curve
    merged_score  = model.predict(data, batch_size=test_batch_size, verbose=0)
    (self.fCurrentAUC, self.fCurrentROCy, self.fCurrentROCx) = GenerateROCCurve(truth, merged_score[:,0])


  ###############################################
  def CreatePlots(self):
    import copy
    """This function saves the plots using the existing results data"""
    AliMLModelResultsBase.CreatePlots(self)

    # Test dataset and evaluate scores
    labels = ['Class 1', 'Class 2']

    ##### Save histograms/plots
    epochs    = range(len(self.fHistoryNumberEvents))

    # AUC
    SavePlot('./Results/{:s}-AUC.png'.format(self.fModelName), 'AUC values', x=epochs, y=(self.fHistoryAUC,), functionlabels=('AUC',), legendloc='lower right', axislabels=('Number of evetns', 'AUC values'))
    # ROC
    SavePlot('./Results/{:s}-ROC.png'.format(self.fModelName), 'ROC curve', x=self.fCurrentROCy, y=(self.fCurrentROCx,self.fCurrentROCy), functionlabels=('(AUC={0:.3f})'.format(self.fCurrentAUC),'Guess ROC'), rangex=(0,1.1), legendloc='lower right', axislabels=('False Positive Rate', 'True Positive Rate') )


###############################################
def SavePlot(fname, title, y, functionlabels, x=[], rangex=(), rangey=(),legendloc='upper right', axislabels=(), logY=False):
  """Simple plot helper function"""
  # Check input data
  for item in x:
    if math.isnan(item):
      print('Plot {:s} was not saved due to invalid values.'.format(fname))
      return
  for i in range(len(y)):
    for item in y[i]:
      if math.isnan(item):
        print('Plot {:s} was not saved due to invalid values.'.format(fname))
        return
  if not len(x):
    x = range(len(y[0]))
  for i in range(len(y)):
    plot.plot(x, y[i], label=functionlabels[i], linewidth=2.0)
  plot.subplots_adjust(left=0.15, right=0.95, top=0.90, bottom=0.15)
  plot.tick_params(labelsize=15, length=10, width=2)
  plot.legend(loc=legendloc, prop={'size': 15})
  if len(rangex) == 2:
    plot.xlim(rangex[0], rangex[1])
  if len(rangey) == 2:
    plot.ylim(rangey[0], rangey[1])
  plot.title(title, fontsize=20, pad=10)
  if len(axislabels) == 2:
    plot.xlabel(axislabels[0], fontsize=15)
    plot.ylabel(axislabels[1], fontsize=15)

#  newfname = fname.replace('.png', '.0.png')
#  i = 1
#  while os.path.isfile(newfname):
#    newfname = fname.replace('.png', '.{:d}.png'.format(i))
#    i += 1
#  fname = newfname

  if logY:
    plot.yscale('log', nonposy='clip')

  plot.savefig("{:s}".format(fname), dpi=128)
  #plot.show()
  plot.clf()


###############################################
def SaveHistogram(fname, title, y, functionlabels, rangex=(0,1), legendloc='upper right', logY=False, nbins=100):
  """Simple histogram helper function"""
  # Check input data
  for i in range(len(y)):
    for item in y[i]:
      if math.isnan(item):
        print('Histogram {:s} was not saved due to invalid values.'.format(fname))
        return


  for i in range(len(y)):
    plot.hist(y[i], label=functionlabels[i], bins=nbins, histtype='step')
  plot.legend(loc=legendloc)
  plot.title(title)
  plot.xlim(rangex[0], rangex[1])
#  plot.yscale('log', nonposy='clip')

  if logY:
    plot.yscale('log', nonposy='clip')


#  newfname = fname.replace('.png', '.0.png')
#  i = 1
#  while os.path.isfile(newfname):
#    newfname = fname.replace('.png', '.{:d}.png'.format(i))
#    i += 1
#  fname = newfname

  plot.savefig("{:s}".format(fname), dpi=128)
  #plot.show()
  plot.clf()


###############################################
def PlotConfusionMatrix(fname, cm, classes, normalize=True, title='Confusion matrix', cmap=None):
  # Code taken from scikit learn documentation
  import itertools

  if cmap is None:
    cmap = plot.cm.Blues

  if normalize:
    cm = cm.astype('float') / cm.sum(axis=0)[:, numpy.newaxis]

  plot.imshow(cm, interpolation='nearest', cmap=cmap)
  plot.title(title)
  plot.colorbar()
  tick_marks = numpy.arange(len(classes))
  plot.xticks(tick_marks, classes, rotation=45)
  plot.yticks(tick_marks, classes)

  fmt = '.2f' if normalize else 'd'
  thresh = cm.max() / 2.
  for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
    plot.text(j, i, format(cm[i, j], fmt), horizontalalignment="center", color="white" if cm[i, j] > thresh else "black")

  plot.ylabel('True label')
  plot.xlabel('Predicted label')

  plot.savefig("{:s}".format(fname), dpi=128)
  plot.clf()



###############################################
def CalculateTaggingRatesBinaryClassifier(model, scores, scores_ref, test_batch_size, eff=1.0, refEff=1.0, verbose=0):
  """Show and return classification efficiencies when demanding different efficiencies in data_ref"""
  scores_ref = copy.deepcopy(scores_ref)
  scores_ref = numpy.sort(scores_ref, axis=0)
  currentThresholds = {}
  for perc in ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9']:
    currentThresholds[perc] = scores_ref[int(float(perc)*(len(scores_ref)))][0]

  # Define the counters and check the scores
  tagged   = {'0.1': 0, '0.2': 0, '0.3': 0, '0.4': 0, '0.5': 0, '0.6': 0, '0.7': 0, '0.8': 0,'0.9': 0}
  for i in range(len(scores)):
    score = scores[i]
    for perc in ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9']:
      if score <= currentThresholds[perc]: # prediction is class 0
        tagged[perc] += 1

  if verbose == 1:
    for perc in ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9']:
      print('At efficiency {:.1f}% (abs: {:2.3f}%)(score < {:E}), tagging rate={:3.4f}%, absolute rate={:3.4f}%'.format(100.*float(perc), 100.*float(perc)*refEff, currentThresholds[perc], 100.*(float(tagged[perc])/float(len(scores))), 100.*(float(tagged[perc])/float(len(scores)))*eff))

  return ([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],[float(tagged[perc])/len(scores) for perc in ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9']])


###############################################
def GenerateROCCurve(truth, score):
  """ROC curve & AUC are generated"""
  from sklearn.metrics import roc_curve, auc
  currentROCy, currentROCx, _ = roc_curve(truth, score)
  currentAUC = auc(currentROCy, currentROCx)
  print('AUC={:f}'.format(currentAUC))
  return (currentAUC, currentROCy, currentROCx)
