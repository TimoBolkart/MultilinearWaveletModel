cd ..\VS2010\LocalMultilinearWavelet

call Release\LocalMultilinearWavelet.exe ..\..\Model\bu3dfe_all_scaled.lmm ..\..\Model\MeanFace_NE00_Scaled.off ..\..\Model\ScaledAlignedLandmarks.txt ..\..\Example\stereo_pointcloud.off ..\..\Example\stereo_pointcloud_landmarks.txt ..\..\Example\stereo_pointcloud_fitting 10.0 0.5 0.0 0.0

pause