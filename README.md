# MSDL_based_artifact_rejection
submitted to Neurocomputing

__\<Despription\>__

These codes use an open-access dataset, <a href="http://www.bbci.de/competition/iv/#datasets" target="_blank">BCI competition IV dataset 2a</a>.
Please send <a href="http://www.bbci.de/competition/iv/#download" target="_blank">e-mail</a> to access the data before using our codes.

This project has two main m.files:
1. Make_ss_data
    - trigg (select biosig4octmat-3.1.0 (<a href="https://sourceforge.net/projects/biosig/files/BioSig%20for%20Octave%20and%20Matlab/" target="_blank">link</a>))
    - sload
    - automatic_wICA
      - runica (<a href="https://sccn.ucsd.edu/eeglab/download.php" target="_blank">
     link</a>)
2. SS_analysis
    - sep_AF
    - sep_EMDandCC
      - emd (<a href="https://jp.mathworks.com/matlabcentral/fileexchange/52502-denoising-signals-using-empirical-mode-decomposition-and-hurst-analysis?focused=5516501&tab=function" target="_blank">link</a>)
    - sep_DOAC
      - optimize_parameters
      - minFunc(<a href="https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html" target="_blank">link</a>)
      - mi(<a href="http://www.cs.man.ac.uk/~pococka4/MIToolbox.html" target="_blank">link</a>)
    - sep_MSDL
      - decompose
        - blockbasedMP
        - multipassMP
        - multistaegMP
        - shiftinvariantMP
      - learning
        - MPSVD
        - multistageWaveformLearning
      - utilities
        - atomsToSources
        - fasterXcorr
        - sourcesToSignalComponents

and six folders:
1. data
    - A01T_label.mat
    - A01T.gdf
   ...
2. minFunc_2012 (after downloading this file, you have to transfer optimize_parameters.m from MSDL_based_artifact_rejection to MSDL_based_artifact_rejection/minFunc_2012)
3. MIToolbox-3.0.0
4. decompose
5. learning
6. utilities.

The last three files (4-6) are downloaded files (<a href="http://cnel.ufl.edu/~ajbrockmeier/eeg/" target="_blank">Algorithms and demos.zip</a>)

__\<Environments\>__
MALBAB R2017a
 1. Signal Processing Toolbox
 2. Statics and Machine Learning Toolbox
 3. Wavelet Toolbox
 4. biosig4octmat-3.1.0
!! if you have optimization toolbox, you can change minFunc to fminunc.
