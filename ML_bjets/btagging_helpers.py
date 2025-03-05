from __future__ import print_function
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LinearRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR
from sklearn.model_selection import cross_val_score

import numpy as np
import pandas as pd
import uproot
import os.path, math
import matplotlib.pyplot as plt
import ROOT
from functools import partial

gTarget = 2


#########################################################
def LoadInputData_MC(recreate, dataname, nJets):

    # Load data from HFD5 store if already created
    data = None
    if os.path.isfile("./Data/{}/dataset.h5".format(dataname)) and not recreate:
        store = pd.HDFStore("./Data/{}/dataset.h5".format(dataname))
        if "MC" in store:
            print("")
            print(
                "### Loading input dataset ({}) from HDF5 store... ####".format("Data")
            )
            data = store["MC"]

    # If data could not be loaded -> create the HDF5 store
    if data is None:
        print("")
        print(
            "#### Creating input dataset for HDF5 store... this will take some time ####"
        )

        # Open the ROOT file and access the TTree
        file = uproot.open("Data/MergedTTree_{}.root".format(dataname))
        tree = file["jetTree"]
        # tree = file["myTree"]

        # Define the branches you want to read
        branches = [
            "jetPt",
            "jetEta",
            "jetPhi",
            "nTracks",
            "nSV",
            "jetMass",
            "jetFlavor",
            "trackPt",
            "trackEta",
            "dotProdTrackJet",
            "dotProdTrackJetOverJet",
            "deltaRJetTrack",
            "signedIP2D",
            "signedIP3D",
            "momFraction",
            "deltaRTrackVertex",
            "svPt",
            "deltaRSVJet",
            "svMass",
            "svfE",
            "svIPxy",
            "svCPA",
            "svChi2",
            "svDispersion",
            "svDecayLength2D",
            "svDecayLength3D"
        ]

        # Convert the TTree to a Pandas DataFrame
        df = tree.arrays(branches, library="pd")

        # Flatten the arrays for track constituents and secondary vertices
        df = df.apply(pd.Series.explode)

        # Filter DataFrames based on jet flavor
        lf_jet_df = df[(df["jetFlavor"] == 0) | (df["jetFlavor"] == 3)]
        c_jet_df = df[df["jetFlavor"] == 1]
        b_jet_df = df[df["jetFlavor"] == 2]

        # Select 1000 rows from each DataFrame
        lf_jet_sample = (
            lf_jet_df.sample(n=int(nJets * 0.8), random_state=1)
            if len(lf_jet_df) >= (nJets * 0.8)
            else lf_jet_df
        )
        c_jet_sample = (
            c_jet_df.sample(n=int(nJets * 0.2), random_state=1)
            if len(c_jet_df) >= (nJets * 0.2)
            else c_jet_df
        )
        b_jet_sample = (
            b_jet_df.sample(n=nJets, random_state=1)
            if len(b_jet_df) >= nJets
            else b_jet_df
        )

        data = pd.concat([lf_jet_sample, c_jet_sample, b_jet_sample], copy=False)

        # Replace infinite values with NaN and then replace NaNs with zeros
        data.replace([np.inf, -np.inf], np.nan, inplace=True)
        data.fillna(0, inplace=True)

        # Save the data to the HDF5 store
        store = pd.HDFStore("./Data/{}/dataset.h5".format(dataname))
        store.put("MC", data)
        # store["MC"] = data

    print("... done")
    print("")
    data = data.fillna(0.0)
    return data


#########################################################
def LoadInputData_Data(recreate):

    # Load data from HFD5 store if already created
    data = None
    if os.path.isfile("./dataset.h5") and not recreate:
        store = pd.HDFStore("./dataset.h5")
        if "Data" in store:
            print("")
            print(
                "### Loading input dataset ({}) from HDF5 store... ####".format("Data")
            )
            data = store["Data"]

    # If data could not be loaded -> create the HDF5 store
    if data is None:
        print("")
        print(
            "#### Creating input dataset for HDF5 store... this will take some time ####"
        )

        # Open the ROOT file and access the TTree
        file = uproot.open("Data/MergedTTree_LHC23d4_80Percent.root")
        tree = file["bjet-tree-merger/myTree"]

        # Define the branches you want to read
        branches = [
            "jetpT",
            "jetEta",
            "jetPhi",
            "nTracks",
            "nSV",
            "jetMass",
            "jetFlavor",
            "trackpT",
            "trackEta",
            "dotProdTrackJet",
            "dotProdTrackJetOverJet",
            "deltaRJetTrack",
            "signedIP2D",
            "signedIP3D",
            "momFraction",
            "deltaRTrackVertex",
            "svpT",
            "deltaRSVJet",
            "svMass",
            "svfE",
            "svIPXY",
            "svCPA",
            "svChi2PCA",
            "svDispersion",
            "svDecayLength2D",
            "svDecayLength3D"
        ]

        # Convert the TTree to a Pandas DataFrame
        mydf = tree.arrays(branches, library="pd")

        # Flatten the arrays for track constituents and secondary vertices
        mydf = mydf.apply(pd.Series.explode)

        # Select 1000 rows from each DataFrame
        data = mydf.sample(n=1000, random_state=1) if len(mydf) >= 1000 else mydf

        # data = mydf[:1000]

        # Save the data to the HDF5 store
        store = pd.HDFStore("./dataset.h5")
        store["Data"] = data

    print("... done")
    print("")
    data = data.fillna(0.0)
    return data


#########################################################
def CheckData(data, name):
    """
    Data exploration
    """

    import seaborn as sns

    print("Checking Data/MC")

    dataC = pd.DataFrame(data)

    # Replace infinite values with NaN and then replace NaNs with zeros
    dataC.replace([np.inf, -np.inf], np.nan, inplace=True)
    dataC.fillna(0, inplace=True)

    dataC.hist(bins=50, figsize=(20, 15))
    plt.savefig("./Results/{}/DataOverview.png".format(name), dpi=120)
    plt.clf()

    # correlation matrix
    sns.heatmap(dataC.corr(), square=True, vmin=-1, vmax=1, cmap="PiYG")
    plt.subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=0.10)
    plt.savefig("./Results/{}/DataCorrelation.png".format(name), dpi=120)
    plt.clf()

    # Pair correlations
    sns.pairplot(
        dataC,
        hue="jetFlavor",
        vars=[
            "signedIP2D[0]",
            "signedIP2D[1]",
            "signedIP2D[2]",
            "signedIP3D[0]",
            "signedIP3D[1]",
            "signedIP3D[2]",
            "svDecayLength2D[0]",
            "svDecayLength2D[1]",
            "svDecayLength2D[2]",
            "svDecayLength3D[0]",
            "svDecayLength3D[1]",
            "svDecayLength3D[2]",
        ],
    )
    plt.subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=0.10)
    plt.savefig("./Results/{}/PairCorrelations.png".format(name), dpi=120)
    plt.clf()


#########################################################
def AddScoresToFiles_MC(models, part):
    """Add scores to file
    Can be used in multiprocessing
    """
    import btagging_keras, btagging_models

    fname = "/Users/hadi/testbjetAnalysis/ML_bjets/TEST_{}.root".format(part)

    # Open the ROOT file and access the TTree
    file = uproot.open(
        "/Users/hadi/testbjetAnalysis/MergedTTree_LHC23d4_80Percent.root"
    )
    tree = file["bjet-tree-merger/myTree"]

    # Define the branches you want to read
    branches = [
        "mJetpT",
        "mJetEta",
        "mJetPhi",
        "mNTracks",
        "mNSV",
        "mJetMass",
        "jetFlavor",
        "mTrackpT",
        "mTrackEta",
        "mDotProdTrackJet",
        "mDotProdTrackJetOverJet",
        "mDeltaRJetTrack",
        "mSignedIP2D",
        "mSignedIP2DSign",
        "mSignedIP3D",
        "mSignedIP3DSign",
        "mMomFraction",
        "mDeltaRTrackVertex",
        "mSVpT",
        "mDeltaRSVJet",
        "mSVMass",
        "mSVfE",
        "mIPXY",
        "mCPA",
        "mChi2PCA",
        "mDispersion",
        "mDecayLength2D",
        "mDecayLength2DError",
        "mDecayLength3D",
        "mDecayLength3DError",
    ]

    # Convert the TTree to a Pandas DataFrame
    df = tree.arrays(branches, library="pd")

    # Flatten the arrays for track constituents and secondary vertices
    df = df.apply(pd.Series.explode)

    # Filter DataFrames based on jet flavor
    data_b = df[df["jetFlavor"] == 2]
    data_c = df[df["jetFlavor"] == 1]
    data_udsg = df[df["jetFlavor"] == 0 or df["jetFlavor"] == 3]

    for dat in [data_b, data_c, data_udsg]:
        dat.drop(["jetFlavor"], inplace=True, axis=1)

    # Load model and apply it
    print("### Adding model predictions ...")
    predictions_b = []
    predictions_c = []
    predictions_udsg = []
    for modelName in models:
        if not "keras" in modelName.lower():
            import joblib

            model = joblib.load(
                "./Models/Model_{}_Target{}.pkl".format(modelName, gTarget)
            )
            predictions_b += [model.predict_proba(data_b)[:, 1]]
            predictions_c += [model.predict_proba(data_c)[:, 1]]
            predictions_udsg += [model.predict_proba(data_udsg)[:, 1]]
        else:
            model = btagging_keras.AliMLKerasModel(2)
            model.LoadModel("Model_{}_Target{}".format(modelName, gTarget))
            data_keras_b = btagging_models.GetDataForKeras(data_b)
            data_keras_c = btagging_models.GetDataForKeras(data_c)
            data_keras_udsg = btagging_models.GetDataForKeras(data_udsg)
            predictions_b += [
                model.fModel.predict(data_keras_b, batch_size=512, verbose=0)[:, 0]
            ]
            predictions_c += [
                model.fModel.predict(data_keras_c, batch_size=512, verbose=0)[:, 0]
            ]
            predictions_udsg += [
                model.fModel.predict(data_keras_udsg, batch_size=512, verbose=0)[:, 0]
            ]
    print("... done")

    # Write results as new Tree to file
    AddBranchesToTreeInFile(
        fname, "Scores_{}".format("bJets"), models, predictions_b, overwrite=True
    )
    AddBranchesToTreeInFile(
        fname, "Scores_{}".format("cJets"), models, predictions_c, overwrite=True
    )
    AddBranchesToTreeInFile(
        fname, "Scores_{}".format("udsgJets"), models, predictions_udsg, overwrite=True
    )
    print("... Added predictions to Trees")


#########################################################
def AddScoresToFiles(models, part):
    """Add scores to file
    Can be used in multiprocessing
    """
    import btagging_keras, btagging_models

    fname = "/home/hadi-hassan/bjetAnalysisML/ML_bjets/test_tree_score_{}.root".format(
        part
    )

    # Open the ROOT file and access the TTree
    file = uproot.open("/home/hadi-hassan/bjetAnalysisML/test_tree.root")
    tree = file["bjet-tree-merger/myTree"]
    # Define the branches you want to read
    branches = [
        "mJetpT",
        "mJetEta",
        "mJetPhi",
        "mNTracks",
        "mNSV",
        "mJetMass",
        "jetFlavor",
        "mTrackpT",
        "mTrackEta",
        "mDotProdTrackJet",
        "mDotProdTrackJetOverJet",
        "mDeltaRJetTrack",
        "mSignedIP2D",
        "mSignedIP2DSign",
        "mSignedIP3D",
        "mSignedIP3DSign",
        "mMomFraction",
        "mDeltaRTrackVertex",
        "mSVpT",
        "mDeltaRSVJet",
        "mSVMass",
        "mSVfE",
        "mIPXY",
        "mCPA",
        "mChi2PCA",
        "mDispersion",
        "mDecayLength2D",
        "mDecayLength2DError",
        "mDecayLength3D",
        "mDecayLength3DError",
    ]
    # Convert the TTree to a Pandas DataFrame
    mydf = tree.arrays(branches, library="pd")
    # Flatten the arrays for track constituents and secondary vertices
    mydf = mydf.apply(pd.Series.explode)
    # Select 1000 rows from each DataFrame
    # data = mydf.sample(n=1000, random_state=1) if len(mydf) >= 1000 else mydf
    data = mydf

    # Load model and apply it
    print("### Adding model predictions ...")
    predictions = []
    for modelName in models:
        if not "keras" in modelName.lower():
            import joblib

            model = joblib.load(
                "./Models/Model_{}_Target{}.pkl".format(modelName, gTarget)
            )
            predictions += [model.predict_proba(data)[:, 1]]
        else:
            model = btagging_keras.AliMLKerasModel(2)
            model.LoadModel("LHC23d4_5_20_Track2Coll/Model_{}".format(modelName))
            data_keras = btagging_models.GetDataForKeras(data)
            predictions += [
                model.fModel.predict(data_keras, batch_size=512, verbose=0)[:, 0]
            ]
    print("... done")

    # Write results as new Tree to file
    AddBranchesToTreeInFile(
        fname, "Scores_{}".format("allJets"), models, predictions, overwrite=True
    )
    print("... Added predictions to Trees")


#########################################################
def AddBranchesToTreeInFile(fname, tname, bNames, bData, offset=0, overwrite=False):
    """Add custom branches to tree in file"""

    ##### Try to read file & tree
    ofile = ROOT.TFile(fname, "RECREATE")
    outTree = ofile.Get(tname)
    if not outTree or overwrite == True:
        outTree = ROOT.TTree(tname, tname)

    if overwrite == False:
        offset = outTree.GetEntries()

    ##### Create branches, if they do not exist
    outBranches = [
        None,
    ] * len(bNames)
    outBuffer = [
        None,
    ] * len(bNames)
    for i, bname in enumerate(bNames):
        outBranches[i] = outTree.GetBranch(bname)
        outBuffer[i] = np.zeros(1, dtype=float)
        if not outBranches[i]:
            outBranches[i] = outTree.Branch(bname, outBuffer[i], "{:s}/D".format(bname))
        else:
            outTree.SetBranchAddress(bname, outBuffer[i])

    ##### Loop through the chain and add raw samples to output tree
    for iData in range(len(bData[0])):

        for iType in range(len(bData)):
            try:
                outBuffer[iType][0] = bData[iType][iData][0]
            except:
                outBuffer[iType][0] = bData[iType][iData]

        outTree.GetEntry(offset + iData)
        outTree.Fill()

    outTree.Write("", ROOT.TObject.kOverwrite)
    ofile.Close()
