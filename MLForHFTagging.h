#ifndef MLRESPONSEHFTAGGING_H_
#define MLRESPONSEHFTAGGING_H_

#include <map>
#include <string>
#include <vector>
#include "TaggingUtilities.h"

// C++ and system includes
#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#else
#include <onnxruntime_cxx_api.h>
#endif

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_BJET(FEATURE)                                 \
  {                                                            \
    #FEATURE, static_cast<uint8_t>(InputFeaturesBTag::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the VECTOR vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_BTAG_FULL(VECTOR, OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesBTag::FEATURE):              \
  {                                                                   \
    VECTOR.emplace_back(OBJECT.GETTER);                               \
    break;                                                            \
  }

// Specific case of CHECK_AND_FILL_VEC_BTAG_FULL(VECTOR, OBJECT, FEATURE, GETTER)
// where FEATURE = GETTER
#define CHECK_AND_FILL_VEC_BTAG(VECTOR, OBJECT, GETTER) \
  case static_cast<uint8_t>(InputFeaturesBTag::GETTER): \
  {                                                     \
    VECTOR.emplace_back(OBJECT.GETTER);                 \
    break;                                              \
  }

struct BJetParams
{
  float jetpT = 0.0;
  float jetEta = 0.0;
  float jetPhi = 0.0;
  int nTracks = -1;
  int nSV = -1;
  float jetMass = 0.0;
};

struct BJetTrackParams
{
  double trackpT = 0.0;
  double trackEta = 0.0;
  double dotProdTrackJet = 0.0;
  double dotProdTrackJetOverJet = 0.0;
  double deltaRJetTrack = -1;
  double signedIP2D =0.0;
  double signedIP2DSign =0.0;
  double signedIP3D = 0.0;
  double signedIP3DSign = 0.0;
  double momFraction = 0.0;
  double deltaRTrackVertex = -1;
  double trackPhi = 0.0;
  double trackCharge = 0.0;
  double trackITSChi2NCl = 0.0;
  double trackTPCChi2NCl = 0.0;
  double trackITSNCls = 0.0;
  double trackTPCNCls = 0.0;
  double trackTPCNCrossedRows = 0.0;
  int trackOrigin = -1;
  int trackVtxIndex = -1;
};

struct BJetSVParams
{
  double svPt = 0.0;
  double deltaRSVJet = -1;
  double svMass = 0.0;
  double svfE = 0.0;
  double svIPxy = 0.0;
  double svCPA = 0.0;
  double svChi2PCA = 0.0;
  double svDispersion = 0.0;
  double svDecayLength2D = 0.0;
  double svDecayLength2DError = 0.0;
  double svDecayLength3D = 0.0;
  double svDecayLength3DError = 0.0;
};

enum class InputFeaturesBTag : uint8_t
{
  jetpT = 0,
  jetEta,
  jetPhi,
  nTracks,
  nSV,
  jetMass,
  trackpT,
  trackEta,
  dotProdTrackJet,
  dotProdTrackJetOverJet,
  deltaRJetTrack,
  signedIP2D,
  signedIP2DSign,
  signedIP3D,
  signedIP3DSign,
  momFraction,
  deltaRTrackVertex,
  trackPhi,
  trackCharge,
  trackITSChi2NCl,
  trackTPCChi2NCl,
  trackITSNCls,
  trackTPCNCls,
  trackTPCNCrossedRows,
  trackOrigin,
  trackVtxIndex,
  svPt,
  deltaRSVJet,
  svMass,
  svfE,
  svIPxy,
  svCPA,
  svChi2PCA,
  svDispersion,
  svDecayLength2D,
  svDecayLength2DError,
  svDecayLength3D,
  svDecayLength3DError,
};

class MLForHFTagging
{
public:
  /// Default constructor
  MLForHFTagging() = default;
  /// Default destructor
  virtual ~MLForHFTagging() = default;

  /// Method to fill the inputs of jet, tracks and secondary vertices
  /// \param jet is the b-jet candidate
  /// \param tracks is the vector of tracks associated to the jet
  /// \param svs is the vector of secondary vertices associated to the jet
  /// \param jetInput the jet input features vector to be filled
  /// \param trackInput the tracks input features vector to be filled
  /// \param svInput the SVs input features vector to be filled
  template <typename T1, typename T2, typename T3>
  void fillInputFeatures(T1 const &jet, std::vector<T2> const &tracks, std::vector<T3> const &svs, std::vector<float> &jetInput, std::vector<float> &trackInput, std::vector<float> &svInput)
  {

    // Jet features
    for (const auto &idx : mCachedIndices)
    {
      switch (idx)
      {
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, jetpT)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, jetEta)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, jetPhi)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, nTracks)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, nSV)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, jetMass)

      default:
        break;
      }
    }

    // Track features
    for (const auto &track : tracks)
    {
      for (const auto &idx : mCachedIndices)
      {
        switch (idx)
        {
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackpT)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackEta)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, dotProdTrackJet)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, dotProdTrackJetOverJet)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, deltaRJetTrack)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, signedIP2D)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, signedIP2DSign)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, signedIP3D)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, signedIP3DSign)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, momFraction)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, deltaRTrackVertex)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackPhi)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackCharge)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackITSChi2NCl)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackTPCChi2NCl)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackITSNCls)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackTPCNCls)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackTPCNCrossedRows)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackOrigin)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackVtxIndex)

        default:
          break;
        }
      }
    }

    // Secondary vertex features
    for (const auto &sv : svs)
    {
      for (const auto &idx : mCachedIndices)
      {
        switch (idx)
        {
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svPt)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, deltaRSVJet)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svMass)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svfE)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svIPxy)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svCPA)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svChi2PCA)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svDispersion)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svDecayLength2D)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svDecayLength2DError)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svDecayLength3D)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svDecayLength3DError)

        default:
          break;
        }
      }
    }
  }

  template <typename T>
  static int replaceNaN(std::vector<T>& vec, T value)
  {
    int numNaN = 0;
    for (auto& el : vec) {
      if (std::isnan(el) || std::isinf(el) || std::abs(el) == 99 || std::abs(el) == 999) {
        el = value;
        ++numNaN;
      }
    }
    return numNaN;
  }


  /// Method to get the input features vector needed for ML inference in a 2D vector
  /// \param jet is the b-jet candidate
  /// \param tracks is the vector of tracks associated to the jet
  /// \param svs is the vector of secondary vertices associated to the jet
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename T3>
  std::vector<std::vector<float>> getInputFeatures2D(T1 const &jet, std::vector<T2> const &tracks, std::vector<T3> const &svs)
  {

    std::vector<float> jetInput;
    std::vector<float> trackInput;
    std::vector<float> svInput;

    fillInputFeatures(jet, tracks, svs, jetInput, trackInput, svInput);

    std::vector<std::vector<float>> inputFeatures;

    replaceNaN(jetInput, 0.f);
    replaceNaN(trackInput, 0.f);
    replaceNaN(svInput, 0.f);
    
    inputFeatures.push_back(jetInput);
    inputFeatures.push_back(trackInput);
    inputFeatures.push_back(svInput);

    return inputFeatures;
  }

  /// Method to get the input features vector needed for ML inference in a 1D vector
  /// \param jet is the b-jet candidate
  /// \param tracks is the vector of tracks associated to the jet
  /// \param svs is the vector of secondary vertices associated to the jet
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename T3>
  std::vector<float> getInputFeatures1D(T1 const &jet, std::vector<T2> const &tracks, std::vector<T3> const &svs)
  {

    std::vector<float> jetInput;
    std::vector<float> trackInput;
    std::vector<float> svInput;

    fillInputFeatures(jet, tracks, svs, jetInput, trackInput, svInput);

    std::vector<float> inputFeatures;

    inputFeatures.insert(inputFeatures.end(), jetInput.begin(), jetInput.end());
    inputFeatures.insert(inputFeatures.end(), trackInput.begin(), trackInput.end());
    inputFeatures.insert(inputFeatures.end(), svInput.begin(), svInput.end());

    replaceNaN(inputFeatures, 0.f);

    return inputFeatures;
  }

  /// Method to translate configurable input-feature strings into integers
  /// \param cfgInputFeatures array of input features names
  void cacheInputFeaturesIndices(std::vector<std::string> const &cfgInputFeatures)
  {
    setAvailableInputFeatures();
    for (const auto &inputFeature : cfgInputFeatures)
    {
      if (mAvailableInputFeatures.count(inputFeature))
      {
        mCachedIndices.emplace_back(mAvailableInputFeatures[inputFeature]);
      }
      else
      {
        throw std::runtime_error(Form("Input feature %s not available! Please check your configurables.", inputFeature.c_str()));
      }
    }
  }

protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    mAvailableInputFeatures = {
        // Jet features
        FILL_MAP_BJET(jetpT),
        FILL_MAP_BJET(jetEta),
        FILL_MAP_BJET(jetPhi),
        FILL_MAP_BJET(nTracks),
        FILL_MAP_BJET(nSV),
        FILL_MAP_BJET(jetMass),

        // Track features
        FILL_MAP_BJET(trackpT),
        FILL_MAP_BJET(trackEta),
        FILL_MAP_BJET(dotProdTrackJet),
        FILL_MAP_BJET(dotProdTrackJetOverJet),
        FILL_MAP_BJET(deltaRJetTrack),
        FILL_MAP_BJET(signedIP2D),
        FILL_MAP_BJET(signedIP2DSign),
        FILL_MAP_BJET(signedIP3D),
        FILL_MAP_BJET(signedIP3DSign),
        FILL_MAP_BJET(momFraction),
        FILL_MAP_BJET(deltaRTrackVertex),
        FILL_MAP_BJET(trackPhi),
        FILL_MAP_BJET(trackCharge),
        FILL_MAP_BJET(trackITSChi2NCl),
        FILL_MAP_BJET(trackTPCChi2NCl),
        FILL_MAP_BJET(trackITSNCls),
        FILL_MAP_BJET(trackTPCNCls),
        FILL_MAP_BJET(trackTPCNCrossedRows),
        FILL_MAP_BJET(trackOrigin),
        FILL_MAP_BJET(trackVtxIndex),

        // Secondary vertex features
        FILL_MAP_BJET(svPt),
        FILL_MAP_BJET(deltaRSVJet),
        FILL_MAP_BJET(svMass),
        FILL_MAP_BJET(svfE),
        FILL_MAP_BJET(svIPxy),
        FILL_MAP_BJET(svCPA),
        FILL_MAP_BJET(svChi2PCA),
        FILL_MAP_BJET(svDispersion),
        FILL_MAP_BJET(svDecayLength2D),
        FILL_MAP_BJET(svDecayLength2DError),
        FILL_MAP_BJET(svDecayLength3D),
        FILL_MAP_BJET(svDecayLength3DError)};
  }

private:
  std::map<std::string, uint8_t> mAvailableInputFeatures; // map of available input features
  std::vector<uint8_t> mCachedIndices;                    // vector of index correspondance between configurables and available input features
};

#undef FILL_MAP_BJET
#undef CHECK_AND_FILL_VEC_BTAG_FULL
#undef CHECK_AND_FILL_VEC_BTAG

#endif // MLRESPONSEHFTAGGING_H_

#include <iostream> // Include iostream for debugging

// Class for loading the ML model with ONNX
template <std::size_t nConst = 10>
class JetModel
{
public:
  JetModel(const std::string &modelPath, DebugLevel debugLevel = DebugLevel::INFO) : mDebugLevel(debugLevel)
  {
    session = new Ort::Session(Ort::Env(ORT_LOGGING_LEVEL_WARNING, "JetModel"), modelPath.c_str(), Ort::SessionOptions{});
  }

  ~JetModel()
  {
    delete session;
  }

  std::vector<float> predict(std::vector<float> &jetFeatures,
                             std::vector<float> &trackFeatures,
                             std::vector<float> &svFeatures)
  {
    // Create input tensors
    Ort::MemoryInfo memInfo = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);

    // Jet features tensor
    std::vector<int64_t> jetInputShape = {1, static_cast<int64_t>(jetFeatures.size())};
    Ort::Value jetInputTensor = Ort::Value::CreateTensor<float>(memInfo, jetFeatures.data(), jetFeatures.size(), jetInputShape.data(), jetInputShape.size());
    log(mDebugLevel, DebugLevel::DEBUG, "Jet input tensor created with shape: [" + std::to_string(jetInputShape[0]) + ", " + std::to_string(jetInputShape[1]) + "]");

    // Track features tensor
    std::vector<int64_t> trackInputShape = {1, static_cast<int64_t>(nConst), static_cast<int64_t>(trackFeatures.size() / nConst)}; // Adjust based on your model
    Ort::Value trackInputTensor = Ort::Value::CreateTensor<float>(memInfo, trackFeatures.data(), trackFeatures.size(), trackInputShape.data(), trackInputShape.size());
    log(mDebugLevel, DebugLevel::DEBUG, "Track input tensor created with shape: [" + std::to_string(trackInputShape[0]) + ", " + std::to_string(trackInputShape[1]) + ", " + std::to_string(trackInputShape[2]) + "]");

    // SV features tensor
    std::vector<int64_t> svInputShape = {1, static_cast<int64_t>(nConst), static_cast<int64_t>(svFeatures.size() / nConst)}; // Adjust based on your model
    Ort::Value svInputTensor = Ort::Value::CreateTensor<float>(memInfo, svFeatures.data(), svFeatures.size(), svInputShape.data(), svInputShape.size());
    log(mDebugLevel, DebugLevel::DEBUG, "SV input tensor created with shape: [" + std::to_string(svInputShape[0]) + ", " + std::to_string(svInputShape[1]) + ", " + std::to_string(svInputShape[2]) + "]");

    // Run the model
    std::vector<const char *> inputNames = {"B0_GV_In", "B1_CNN1D", "B2_CNN1D"}; // Adjust based on your model's input names
    std::vector<Ort::Value> inputTensors;                                      // Declare the vector to hold input tensors
    inputTensors.push_back(std::move(jetInputTensor));                         // Move the tensor into the vector
    inputTensors.push_back(std::move(trackInputTensor));                       // Move the tensor into the vector
    inputTensors.push_back(std::move(svInputTensor));                          // Move the tensor into the vector

    std::vector<const char *> outputNames = {"dense"}; // Adjust based on your model's output name
    log(mDebugLevel, DebugLevel::DEBUG, "Running the model...");

    try
    {
      Ort::RunOptions runOptions;
      auto outputTensors = session->Run(runOptions, inputNames.data(), inputTensors.data(), inputTensors.size(), outputNames.data(), outputNames.size());
      log(mDebugLevel, DebugLevel::DEBUG, "Model run completed.");

      // Process output
      float *outputData = outputTensors.back().GetTensorMutableData<float>();
      std::vector<float> results(outputData, outputData + 3);
      log(mDebugLevel, DebugLevel::DEBUG, "Output processed. Result: " + std::to_string(results.front()));

      return results;
    }
    catch (const Ort::Exception &e)
    {
      std::cerr << "Error during model run: " << e.what() << std::endl;
      throw;
    }
  }

private:
  DebugLevel mDebugLevel = DebugLevel::INFO;
  Ort::Session *session = nullptr;
};
