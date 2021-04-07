### Parameters for FISH_pipeline.py script

#########################################################################################################
#########################################################################################################
import os
import sys

#########################################################################################################
#########################################################################################################
##### input/output files/Users/u.gowthaman/curie/u.gowthaman/analysis/smFISH_analysis/20190927_SMY3014_2_Glu_PP7-Cy3_MS2-Cy5_1/20190927_2-Pos_000_000_mask_cell+nuc+spots.tif


folderIn = "/DATA/lenstra_lab/u.gowthaman/data/smFISH_data/"

fileListIn = [
#              "20190922/SMY3014_gal_1",
#              "20190927/SMY3014_1_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20190927/SMY3014_2_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20190928/YTL1022_1_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20190928/YTL1022_2_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20190928/YTL1023_1_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20190928/YTL1023_2_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20190928/YTL1072_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20191006/SMY3014_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20191006/YTL1017_2_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20191006/YTL1018_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20191006/YTL1020_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20191006/YTL1021_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20191022/SMY3014_1_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20191022/YTL1018_1_Glu_PP7-Cy3_MS2-Cy5_1_1",
#              "20191022/YTL1018_2_Glu_PP7-Cy3_MS2-Cy5_1_1",
#              "20191022/YTL1019_1_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20191022/YTL1019_2_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20191022/YTL1021_1_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20191022/YTL1023_2_Glu_PP7-Cy3_MS2-Cy5_1",
#              "20191022/YTL1023_3_Glu_PP7-Cy3_MS2-Cy5_1"
#              "20191117/SMY3014_1_GAL_PP7-Cy3_MS2-Cy5_1",
#              "20191117/SMY3014_2_GAL_PP7-Cy3_MS2-Cy5_1",
#              "20191117/YTL1017_1_GLU_PP7-Cy3_MS2-Cy5_1",
#              "20191117/YTL1017_2_GLU_PP7-Cy3_MS2-Cy5_1",
#              "20191117/YTL1020_GLU_PP7-Cy3_MS2-Cy5_1",
#              "20191117/YTL1021_GLU_PP7-Cy3_MS2-Cy5_1",
#              "20191117/YTL1022_GLU_PP7-Cy3_MS2-Cy5_1",
#              "20191117/YTL1072_1_GLU_PP7-Cy3_MS2-Cy5_1",
#              "20191117/YTL1072_2_GLU_PP7-Cy3_MS2-Cy5_1",
#              "20191117/YTL1072_1_GAL_PP7-Cy3_MS2-Cy5_1",
#              "20191117/YTL1072_2_GAL_PP7-Cy3_MS2-Cy5_1"
              "20191127/YTL1017_1_Glu_PP7-Cy3_MS2-Cy5_1",
              "20191127/YTL1017_2_Glu_PP7-Cy3_MS2-Cy5_1",
              "20191127/YTL1019_1_Glu_PP7-Cy3_MS2-Cy5_1",
              "20191127/YTL1020_Glu_PP7-Cy3_MS2-Cy5_1",
              "20191127/YTL1021_Glu_PP7-Cy3_MS2-Cy5_1"
              ]

outputfolder = "/DATA/lenstra_lab/u.gowthaman/analysis/smFISH_analysis/"
PipelineCellProfiler = ""
analysisfolder = "/DATA/lenstra_lab/u.gowthaman/analysis/smFISH_analysis/"

######################################
##### pipeline #######################
######################################
processAll = 1
nrFilesToAnalyze = 0
channelsToAnalyze = [0,1] ### list of channels to analyze in: 0 = Cy3, 1 = Cy5.

MaxProj = 0
RunCellprofiler = 0
RunOptimizeThresh = 0
RunFindCells = 0
RunLocalize = 0
MakeMaskImage = 0
CalculateTranscriptionSites = 1
MakeHistograms = 1
CalcSingleCellCorr = 0
Violins = 0 #making violin plots only possible if also calculating histograms
CompareNormValues = 0
CombineReps = 0 #combine replicates and make histograms and violin plots (optional) for those
UsePrecalcDistr = 0 #Set to 1 if you want to use precalculated datasets (if they exist)


######################################
##### Parameters Max projection ######
######################################
zSlices = 13        # nr z slices
nChannels = 3       # nr channels including DAPI channel
totalImg = zSlices * nChannels

######################################
##### Parameters FindCells ######
######################################
channelCellMask = "Cy3" #choose Cy3 or Cy5 for cell segmentation
ccdist = 25 #(approximate) minimum distance between the centers of cells, good values: yeast: 25, mamalian: 150
threshold = None #pixel value used to threshold the image, default:Otsu
thresholdnuc = None #pixel value used to threshold the image with nuclei, default:Otsu
removeborders = True #remove any cells (and their nuclei) touching any borders (case with nan's unhandled atm), default: True

######################################
##### Parameters Localize ############
######################################
localizeDim = "3D" # either "2D" or "3D". For 3D, image with 0.3 um z step and take more frames above and below cell for fitting.
threshold_cy3=35 # threshold for Cy3 channel
threshold_cy5=18 # threshold for Cy5 channel
psfPx=2.18 # PSF width in pixels. Can be calculated by: (wavelength emitted light * magnification) / (2 * NA objective * camera pixel size). NA objective 40x = 1.4. camera pixel size = 6500 nm. With 40x objective with optovar 1.25, psfPx = 1.70 (cy3 = 1.57 /cy5 = 1.84), for 40x objective with optovar 1.6, pxfPx = 2.18 (cy3 = 2.0 /cy5 = 2.36)
psfPxZ = 0.93 # PSF in z in pixels. Can be calculated by: (wavelength * RI of mounting media) / ((NA objective)^2 * zstep). RI  = 1.47 for prolong gold. NA objective 40x = 1.4. For 40x objective, 500 nm zstep, psfZPx = 0.93 (cy3 = 0.855, cy5 = 1.005)
maxDist=3.; # Maximal distance tolerated between guess and fit (in PSF width)
minSeparation=.1; # Minimal distance tolerated between two spots (in PSF width)
winSize = 5 # size of the window used for PSF fitting
winSizeZ = 3 # size of the window used for PSF fitting in z


######################################
##### Parameters Calculate TS ########
######################################
MaskTS = "nucleus" # mask used for defining the TS (brightest spot). choose between "nucleus" or "cell"
AreaNorm = "cyto" # area used for the normalizing spots. Choose between "cyto" (in cell but not in nucleus), "nucleus" or "cell" (entire cell including nucleus).

ExclEmpty = 0 #excludes cells that are completely empty, no RNAs in either nucleus or cytoplasm. Uses channel that is analyzed to determine if cell is empty. If you want to use another channel, set ExlcEmpty 0 and set parameter useControlForEmptyCells to 1.
useControlForEmptyCells = 0 #uses control hybridization (can be other channel) to determine whether cell is empty (no RNA in cell or nucleus). Removes these cells from analysis in all colors. Set control channel with parameter controlColor. Can only be used if ExclEmpty = 0
controlColor = 0 #0 = Cy3, 1 = Cy5. Set which probe set is control. Is used if useControlForEmptyCells is 1.
useFit = 0 # if 1: normalize TS distrubtion by intensity from fit of the cytoplasmic intensity distribution with a gaussian model. If 0: normalize TS distribution by median intensity cytoplasmic spots

includeMultTS = 0 # 0 for one TS per nucleus, 1 for multiple TS per nucleus using threshold (specify thresholdTS), 2 for multiple TS per nucleus by fixed amount (specify nrTS2Include).
thresholdTS = 1.5 # threshold used to determine if nuclear spot is TS. The hreshold indicates how many nascent RNAs need to present to qualify as TS (it is is multiplied by the normalization factor of the cytoplasmic spot intensity). Note, the results do not contain the cells without TS.
nrTS2Include = 2 # number of TS to include in the analysis. It will select the brightest spots as TS. Note, if less TS are found in a nucleus, the script will add a zero (so each cell will have same nr of TS).


######################################
##### Parameters histogram #########
######################################

exclBins=0  # nr bins to exclude for normalizing and plotting (and fitting) TS intensity histograms (1 means that the first (zero) bin will be exluded, 2 means that the first 2 (0 and 1) bins will be excluded, etc).
exclBinsCount = 0 # nr bins to exclude for plotting and calculating statistics on the cell and nuclear count histograms (1 means that the first (zero) bin will be exluded, 2 means that the first 2 (0 and 1) bins will be excluded, etc). This value is not used for the TS count distribution.

FitBurst = 1 #fitting models is only possible if also calculating histogram
FitPoisson = 1 #fitting models is only possible if also calculating histogram
FitBurstCombi = 0
FitPoissonCombi = 0
CalcFracOff = 0 # calculates the fraction of cells that are transcriptionally silent; only relevant if exclBins = 0, exclBinsCount = 0 and ExclEmpty = 1

### nr bins for cytoplasmic and nuclear/TS histogram plots
nrbins_cy3 = 25 # for cellular and nuclear count hist
nrbins_cy5 = 25 # for cellular and nuclear count hist
nrbins_cy3_ts = 10 # for TS intensity hist
nrbins_cy5_ts = 10 # for TS intensity hist

### set scaling of frequency (y) axis in histograms. If "None", scaling will be done automatically
freqAxisCell_cy3 = 1 # for cellular count hist
freqAxisCell_cy5 = 1 # for cellular count hist
freqAxisNuc_cy3 = 1 # for nuclear count hist
freqAxisNuc_cy5 = 1 # for nuclear count hist
freqAxisTS_cy3 =  1 # for TS intensity hist
freqAxisTS_cy5 =  1 # for TS intensity hist

######################################
##### Parameters CalcSingleCellCorr ##
######################################
filterOnDistance = 0 # set to 1 if you want to filter on colocalized Cy3 - Cy5 spots, set threshold in
distThresh = 4 # threshold for determining colocalization in pixels
onoffThreshCy3 = 6 ### threshold to determin if cells are on or off in Cy3
onoffThreshCy5 = 4 ### threshold to determin if cells are on or off in Cy5


######################################
##### Parameters for Violin plots  ###
######################################
pvalsViolins = [[1,2],[3,4,5,6]] #indicates between which datasets you want to calculate pvalues to be plotted in Violin plots. Starts numbering at 1. Leave empty ([[]]) if you don't want any pvalues calculated.

### set scaling of y-axis in Violin plots. If "None", scaling will be done automatically
yAxisViolinCell_cy3 = 10 # for cellular violin plot Cy3
yAxisViolinCell_cy5 = 10 # for cellular violin plot Cy5
yAxisViolinNuc_cy3 = 8 # for nuclear violin plot Cy3
yAxisViolinNuc_cy5 = 8 # for nuclear violin plot Cy5
yAxisViolinTS_cy3 = 6 # for TS violin plot Cy3
yAxisViolinTS_cy5 = 6 # for TS violin plot Cy5

######################################
##### Parameteres combine reps ######
######################################
repsToCombine = [] #List how you want your replicates to be combined. Format as (for example) [[1,2,3],[4],[5,6]] etc. Make sure all dataset numbers are in this list; numbering starts at 1.
