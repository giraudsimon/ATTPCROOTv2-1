Set(LIBRARY_NAME AtTools)
Set(LINKDEF ${LIBRARY_NAME}LinkDef.h)

set(SRCS
  AtParsers.cxx
  AtELossManager.cxx
  AtTrackTransformer.cxx
  
  AtCSVReader.cxx
  AtEulerTransformation.cxx
  AtSpaceChargeModel.cxx
  AtLineChargeModel.cxx
  AtRadialChargeModel.cxx
  AtEDistortionModel.cxx
  AtVirtualTerminal.cxx
  AtKinematics.cxx
  
  AtFormat.cxx
  AtHitSampling/AtSample.cxx
  AtHitSampling/AtSampleMethods.cxx
  AtHitSampling/AtIndependentSample.cxx
  AtHitSampling/AtChargeWeighted.cxx
  AtHitSampling/AtUniform.cxx
  AtHitSampling/AtSampleFromReference.cxx
  AtHitSampling/AtGaussian.cxx
  AtHitSampling/AtWeightedGaussian.cxx
  AtHitSampling/AtWeightedY.cxx

  ElectronicResponse/AtElectronicResponse.cxx
  ElectronicResponse/AtROOTresponse.cxx
  )

Set(DEPENDENCIES
  ROOT::XMLParser
  ROOT::Core
  
  FairRoot::FairTools
  ATTPCROOT::AtData
  ATTPCROOT::AtParameter
)

Set(INCLUDE_DIR
  ${CMAKE_SOURCE_DIR}/AtTools/AtHitSampling
  ${CMAKE_SOURCE_DIR}/AtTools/ElectronicResponse
  )

generate_target_and_root_library(${LIBRARY_NAME}
  LINKDEF ${LINKDEF}
  SRCS ${SRCS}
  DEPS_PUBLIC ${DEPENDENCIES}
  INCLUDE_DIR ${INCLUDE_DIR}
  )

