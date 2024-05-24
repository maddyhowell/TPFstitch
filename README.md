# TPFstitch

By inputting the FITs files of between 2-9 K2 target pixel files (TPFs) which share a border, TPFstitch produces a single FITs file containing the stitched TPFs, which we call the superstamp. The world coordinate system is maintained within the final product, i.e. the RAs and DECs of targets can be identified in the resultant superstamp. Additionally, the superstamp FITs file can be read straight into asteroseismic data analysis tools, such as LightKurve. 

Currently, TPFstitch is only designed for the NGC 6333 and NGC 6273 clusters, but it could be easily adapted for any K2 object with several TPFs; e.g. other clusters, asteroids/comets, microlensing events, galaxies, etc. There is also the potential for TPFstitch to be applied to other photometric surveys, such as TESS.

![tpfstitch_eg](https://github.com/maddyhowell/TPFstitch/assets/53502531/979a3a67-f228-4da4-8504-dd0e72c1b97f)
