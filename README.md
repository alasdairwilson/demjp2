Code base to generate helioviewer compatibale jpeg2000 files from AIA derived DEM data.

requires: demreg, SSW

hv_dem2jp2		: can be used to produce a single jp2000 file when given an input of an array of DEM values along with a configuration structure which, as well as the details that are required for the header, provides the datamin and max scaling.
dem_details		: is necessary for the use of hv_dem2jp2, this creates the data structure. For an example see batch_dem_jp2
dem_compare_aia		: calculate expected DN values from a supplied DEM.
hv_xml_compliance	: this is taken from the helioviewer/jp2gen repository to make the xml header creation process easier.

batch_dem_jp2		: Best use is to use batch_dem_jp2 to calculate dem maps for individual or runs of times.

----------------------------------------------------------
Input
----------------------------------------------------------
required inputs	:
start_time	: the initial time, in any anytim compatible format.
cadence		: the cadence in seconds of the observations, must be a multiple of 120 if using synoptic data or 12 if using full AIA
n_obs		: number of observations, 1 for a single set of jpegs
fits_base_dir	: directory for saving the fits files, these are not deleted (though the uncompressed versions are through aia_prep). Within this folder the fits files are saved in ./yyyy/mm/dd/<synoptic file name>
jp2_base_dir	: directory for output jp2k files, within this directory the jpegs are stored in ./yyyy/mm/dd/, these directories are created automatically as long as the idl session has the correct permissions.

optional inputs:
get_fits	: if the fits files are already present they will not be redownloaded but just in case there is a get_fits keyword, it will only attempt to download the fits files if this flag is set. default=False
min_snr		: the minimum snr that is required for the DEM to be computed for this pixel. i.e. if the resultant DEM is expected to have poor SNR then it will not be inverted. default=3
sat_lvl		: by default uses AIA saturation level but can be specified by username, one value used for all aia channels.
serr_per	: systematic error assigned to AIA data as a percentage. Despite AIA errors being computed as per the instrument paper the recovered DEM will be fairly poor unless the data is assigned an additional systematic error. default=0 recommended=10
gauss_stdev	: controls the normalisation function (deprecated)
gauss_mean	: controls the normalisation function (deprecated) these have been removed in favour of hardcoded normalisation function in order to facilitate the transition to a 6&7 channel response.
----------------------------------------------------------
Output
----------------------------------------------------------
Outputs 7 jp2000 files with bytescaled values for the Differential Emissiom Measure in evenly spaced bins of width 0.2 between logT 5.7 K and logT 7.1 K. The details of the bytescaling and units of these jp2k files can be found in the xml header. 

The header of the jp2000 files contain all the relevant information about the DEM along with a copy of the fits header from the AIA fits file. In particular from the 94A channel of AIA. The duplication of this information allows the use of the aia header for operational paramters such as exposure times, observation time and the ACS control system mode for systematically ignoring calibration maneouvres.

----------------------------------------------------------
Current Issues
----------------------------------------------------------
In v1.0/v1.01 the hottest areas of the Sun return poor DEMs, either they are wiped due to negative values or they poorly recover the DN values. The issue is that it is unlikely that a set of parameters exist that allow the inversion to simultaneously recover cool and hot areas due to the multitemp responses of 94 in particular (but 131 is a culprit too).

Skipped data is just skipped. Currently the aia synoptic seems to miss some channels during 

----------------------------------------------------------
Change log: (version of the code used can be found in the jp2000 xml header.)
----------------------------------------------------------
v1.0
----------------------------------------------------------
initial github push.

v1.01
----------------------------------------------------------
Negative DEM values correctly skipped. These will appear as 0 in the jp2k bytescale.

v1.1
----------------------------------------------------------
First substantial update, this switches from a universal 6 channel dem response to a 6&7 channel. The additional channel is the splitting into warm and hot components for the 94A channel.
The hope is that the hottest areas can be recovered more accurately.