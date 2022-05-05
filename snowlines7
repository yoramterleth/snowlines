// snowlines version 6

// This script determines the snow cover extent on glaciers in the GLIMS database 
// the topographic correction relies on the USGS 3DEP DEM, so glaciers outside the US 
// require the users to input different elevation data. 

// the method followed is largely based on:

// Rastner, P., Prinz, R., Notarnicola, C., Nicholson, L., Sailer, R.,
// Schwaizer, G. and Paul, F., 2019. 
// On the automated mapping of snow cover on glaciers and calculation 
// of snow line altitudes from multi-temporal landsat data. 
// Remote Sensing, 11(12), p.1410.



// load data //////////////////////////////////////////////////////////////////////////////

// THE GLIMS POLYGON: changing the glacier name should allow the script to run for any glacier in the GlIMS database.
var GLIMS = ee.FeatureCollection('GLIMS/current').filter(ee.Filter.bounds(AOI)).filter('glac_name=="Turner Glacier"');  
// paint polygon to image
var imageGLIMS = ee.Image().float().paint(GLIMS, 'area');

// THE LANDSAT 8 DATA 

// select only days in summer, since ther is no snowlines in winter 
var startDay = 172 // june 1st 
var endDay = 274 // oct 1st 

var startyear = '2020' 
var endyear = '2021' 

var LS8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filter(ee.Filter.bounds(GLIMS))
    .filterDate(startyear,endyear)  
  .filter(ee.Filter.dayOfYear(startDay, endDay));
   

// cloud masking: here we make a new collection with cloud masks to apply later 
var maskL8sr = function(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0); 
  // This sees if there was a saturation flag for each band (0 = no saturation, 1 = saturated data)
  var opticalBands = ee.Image.constant(1).clip(GLIMS) 
  var thermalBands = ee.Image.constant(1).clip(GLIMS) 
  // assign 1 value to all cells flagged as affected by clouds, producing cloud cover map. 
  
    // Replace the original band data with the scaled data and apply the masks.
    return image.addBands(opticalBands, null, true)
      .addBands(thermalBands, null, true)
      .updateMask(qaMask)
      .updateMask(saturationMask);
}

// apply the mask function 
var cloud_cover = LS8.map(maskL8sr).select('SR_B5')

// define radiometric corection function
var rad_corr = function(image){
  return image.multiply(0.0000275).add(-0.2);
};

// clip the collection to area of interest and band of interest 
// var LS_collection = LS8.select('SR_B5').map(rad_corr); //

print('Image collection:', LS8)

// return the acquisition dates of all the images
var dates = LS8
    .reduceColumns(ee.Reducer.toList(), ["system:time_start"])
    .get('list'); 
print('dates',dates); 


// DIGITAL ELEVATION MODEL 

// elevation data from USGS 
var DEM = ee.Image('USGS/3DEP/10m') ; // .filter(ee.Filter.bounds(AOI)) ; 

// calculate slope (in radians)
var slope = ee.Terrain.slope(DEM.select('elevation')).divide(180).multiply(Math.PI) ; //.median()); 
var aspect = ee.Terrain.aspect(DEM.select('elevation')).divide(180).multiply(Math.PI)  ; 

// regridding at Landsat resolution:

// Get information about the LANDSAT projection.
var LSProjection = LS8.median().projection();
print('Landsat projection:', LSProjection);

// depending on required output format, the following steps can regrid the LS data 
// to get more workable snow maps in large AOIs.

var slope_regrid = slope
    // // Force the next reprojection to aggregate instead of resampling.
    // .reduceResolution({
    //   reducer: ee.Reducer.mean(),
    //   maxPixels: 196530// 1024
    // })
    // // Request the data at the desired scale while keeping LS projection.
    // .reproject({
    //   crs: LSProjection
    //   scale: 100
    // });
    
var aspect_regrid = aspect
    // // Force the next reprojection to aggregate instead of resampling. (unnecesarry here)
    // .reduceResolution({
    //   reducer: ee.Reducer.mean(),
    //   maxPixels: 196530// 1024
    // })
    // // Request the data at the desired scale while keeping LS projection.
    // .reproject({
    //   crs: LSProjection
    //   scale:100
    // });
    
///////////////// FUNCTION DEFINITIONS /////////////////////////////////////////////////

/// TOPOGRAPHIC CORRECTION ///

var topo_corr = function(image){
// get solar azimuth and zenith angles from landsat metadata
 var solar_azimuth = ee.Number(image.get('SUN_AZIMUTH')).multiply(Math.PI/180);  
 var solar_zenith = ee.Number(90).subtract((image.get('SUN_ELEVATION'))).multiply(Math.PI/180);
 
 // compute the cos i value for image following Wu et al. 2016 
 var cosiI = ee.Image.constant(solar_zenith).cos().multiply(slope_regrid.cos())
 var cosiII =  ee.Image.constant(solar_zenith).sin().multiply(slope_regrid.sin())
 var cosiIII = (ee.Image.constant(solar_azimuth).subtract(aspect_regrid)).cos()
 var cosiIV = cosiII.multiply(cosiIII)
 var cosi = cosiI.add(cosiIV)
 
 // the Minnaert constant (can be tweaked for different topo corrections)
 var k = ee.Image.constant(1)
 
 // now correct the radiance 
 var topo_corr =  image.multiply(slope_regrid.cos().divide(cosi)).multiply(k)
 
 return topo_corr }
 
 
/// OTSU /////

// define function using Otsu's method of identifying optimal separation threshold 
// this work is from Nicholas CLinton: 
// https://medium.com/google-earth/otsus-method-for-image-segmentation-f5c48f405e#:~:text=Otsu's%20method%20is%20a%20means,2016).
var otsu = function(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  
  var indices = ee.List.sequence(1, size);
  
  // Compute between sum of squares, where each mean partitions the data.
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)));
  });
  
  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
};


// Definition of functions to loop over images in collection //////////////

var snowmapper = function(image){

// first we make an image histogram 
  var histogram = image.reduceRegion({
     reducer: ee.Reducer.histogram({
       maxBuckets: 256 , 
     }), 
     geometry: GLIMS, 
     scale: 30, 
     maxPixels: 10e9, 
   });
  
  // use Otsu's method to find the optimal threshold between snow and ice 
  var threshold = otsu(histogram.get('SR_B5'));

  // use this threshold to determine the snowmap
  var snow_map = image.updateMask(image.gt(threshold)).updateMask(image.lt(1))
        // set all values to 1 (instead of keeping the original values)
        .gt(threshold).selfMask()
        // rename the band
        .rename('snowmap') ; 
  
return snow_map ; 

}; 

//////////////////////////////////////////////////////////////////////////////////////////////////
// Apply functions to collection

// topographical corrections 
var LS_topo_corrected = LS8.map(topo_corr)

// clip the collection to area of interest and band of interest 
var LS_collection = LS_topo_corrected.select('SR_B5') //.map(rad_corr); 

// compute the snowmaps 
var snowmaps = LS_collection.map(snowmapper); 

print('snowmap:',snowmaps); 

//////////////////////////////////////////////////////////////////////////////////////////////////
// export variables to drive 

// export acquisition dates  
Export.table.toDrive({
collection: dates,
description: 'acquisition_dates',
fileFormat: 'CSV'
}); 


// export actual imagery: turn on only when necessary, significant slowdown

// using tools obtained here:  https://github.com/fitoprincipe/geetools-code-editor  

// load external export function 
 // var batch = require('users/fitoprincipe/geetools:batch'); 

// batch.Download.ImageCollection.toDrive(snowmaps, 'snowmaps_final', 
//   {scale: 30, 
//   region:GLIMS, 
//   type: 'float'}); 

// batch.Download.ImageCollection.toDrive(cloud_cover, 'cloudmaps_final', 
//   {scale: 30, 
//   region:GLIMS, 
//   type: 'float'}); 



