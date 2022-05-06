# snowlines
### Determines snow cover on glaciers from Landsat imagery. 


These scripts determine the snow cover extent on glaciers in the GLIMS database and are designed to be run in the google earth engine code editor. 
The method followed is largely based on:

Rastner, P., Prinz, R., Notarnicola, C., Nicholson, L., Sailer, R., Schwaizer, G. and Paul, F., 2019. 
On the automated mapping of snow cover on glaciers and calculation of snow line altitudes from multi-temporal landsat data. 
Remote Sensing, 11(12), p.1410.

- snowlines7.js determines snow cover extent on a glaicer from individual LS images.
- annual_AAR.js is based on snowlines7, but combines imagery from a temproal range into yearly minimum values of surface relfectance. This allows the automated determination of yearly Accumulation Area Ratios.  
