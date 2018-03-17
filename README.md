# MSDL_based_artifact_rejection
submitted to Neurocomputing<br>

__\<Despription\>__ <br>
These codes use an open-access dataset, <a href="http://www.bbci.de/competition/iv/#datasets" target="_blank">BCI competition IV dataset 2a</a>.<br>
Please send <a href="http://www.bbci.de/competition/iv/#download" target="_blank">e-mail</a> to access the data before using our codes.<br>

This project has two main m.files:<br>
1. Make_ss_data<br>
    - trigg (select biosig4octmat-3.1.0 (<a href="https://sourceforge.net/projects/biosig/files/BioSig%20for%20Octave%20and%20Matlab/" target="_blank">link</a>))<br>
    - sload<br>
    - automatic_wICA<br>
      - runica (<a href="https://sccn.ucsd.edu/eeglab/download.php" target="_blank">
     link</a>)<br>
2. SS_analysis<br>
    - sep_AF
    - sep_EMDandCC<br>
      - emd (<a href="https://jp.mathworks.com/matlabcentral/fileexchange/52502-denoising-signals-using-empirical-mode-decomposition-and-hurst-analysis?focused=5516501&tab=function" target="_blank">link</a>)<br>
    - sep_DOAC<br>
      - optimize_parameters<br>
      - minFunc(<a href="https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html" target="_blank">link</a>)<br>
      - mi(<a href="http://www.cs.man.ac.uk/~pococka4/MIToolbox.html" target="_blank">link</a>)<br>
    - sep_MSDL<br>
      - decompose<br>
        - blockbasedMP<br>
        - multipassMP<br>
        - multistaegMP<br>
        - shiftinvariantMP<br>
      - learning<br>
        - MPSVD<br>
        - multistageWaveformLearning<br>
      - utilities<br>
        - atomsToSources<br>
        - fasterXcorr<br>
        - sourcesToSignalComponents<br>

and six folders:<br>
1. data<br>
   - A01T.mat<br>
   - A01T.gdf<br>
   ...<br>
2. minFunc_2012 (after downloading this folder, you have to replace optimize_parameters.m which is including in this master file)<br>
3. MIToolbox-3.0.0<br>
4. decompose<br>
5. learning<br>
6. utilities.<br>

The last three folders (4-6) are downloaded files (<a href="http://cnel.ufl.edu/~ajbrockmeier/eeg/" target="_blank">Algorithms and demos.zip</a>)<br>

__\<Environments\>__ <br>
MALBAB R2017a<br>
 1. Signal Processing Toolbox<br>
 2. Statics and Machine Learning Toolbox<br>
 3. Wavelet Toolbox<br>
 4. biosig4octmat-3.1.0<br>
!! if you have optimization toolbox, you can change minFunc to fminunc.
