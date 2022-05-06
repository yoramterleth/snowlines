// annual_AAR //  

// This script is based on snowlines, and methodology in Rastner et al.(2016) . Rather than producing maps of snowcover form individual landsat images, 
// it bases its snowmaps on temporal stacks of LS imagery, yielding yealry minimum reflectance values. This produces an estimate of the minimum yealry snowcover, 
// making an estimation of the accumulation area ratio (AAR) possible. 

// The topographic correction here relies on Arctic DEM: can be run for any glacier in the arctic with topo correct, 
// or any glacier worldwide without topo correction. When using level 2 products as input imagery the script is thought to be widely applicable. 

// Has been tested with LS8 (Band 5) and LS5 (Band 4) data. 

// Yoram Terleth 2022. 

// load data //////////////////////////////////////////////////////////////////////////////

// THE GLIMS POLYGON: changing the glacier name should allow the script to run for any glacier in the GlIMS database.
var GLIMS = ee.FeatureCollection('GLIMS/current').filter(ee.Filter.bounds(AOItarfala)).filter('glac_name=="Storglaciaren"');
print('GLIMS:', GLIMS)
// paint polygon to image
var imageGLIMS = ee.Image().float().paint(GLIMS, 'area');

// THE LANDSAT 8 DATA 

// select only days in summer, since ther is no snowlines in winter 
var startDay = 180 // june 1st 
var endDay = 271 // oct 1st 

var startyear = "2013" //ee.Date('2018')
var endyear = "2021" //ee.Date('2021')

var startyearDate = ee.Date(startyear)
var endyearDate = ee.Date(endyear)

var LS8_all = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') // ('LANDSAT/LT05/C02/T1_L2') //  
    .filter(ee.Filter.bounds(GLIMS))
    .filterDate(startyear,endyear)  
    .filter(ee.Filter.dayOfYear(startDay, endDay));
    

////////////////////////////////////////////////////////////////////////////////////

// Get yearly minimum values, equivalent to lowest snow extent. 
var step = 1 
var numberOfDays = endyearDate.difference(startyearDate, 'years')

var LS8 = ee.ImageCollection(
  ee.List.sequence(0, numberOfDays.subtract(step))
    .map(function (dayOffset) {
      var start = startyearDate.advance(dayOffset, 'years')
      var end = start.advance(step, 'years')  // temporal resolution needs adjusting 
      return LS8_all.filterDate(start,end).min()
     
    })); 
        


//var LS8 = LS8_all ; 
print(LS8) ; 
/////////////////////////////////////////////////////////////////////
   

// cloud masking from assignment only here we make a new collection with cloud masks to apply later 
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
// var LS_collection = LS8.map(rad_corr); //

print('Image collection:', LS8)

// return the acquisition dates of all the images
var dates = LS8
    .reduceColumns(ee.Reducer.toList(), ["system:time_start"])
    .get('list'); 
print('dates',dates); 


// DIGITAL ELEVATION MODEL 

// elevation data from USGS 
var DEM = ee.ImageCollection('UMN/PGC/ArcticDEM/V3/2m').filter(ee.Filter.bounds(AOIall)).select('elevation').mean()
// calculate slope (in radians)
var slope = ee.Terrain.slope(DEM).divide(180).multiply(Math.PI) ; //.median()); 
var aspect = ee.Terrain.aspect(DEM).divide(180).multiply(Math.PI)  ; 

// cehck whether arcti dem exists here: 
print('Slope:',slope)

// regridding at Landsat resolution:

// Get information about the LANDSAT projection.
var LSProjection = LS8.median().projection();
print('Landsat projection:', LSProjection);

// Get the forest cover data at MODIS scale and projection.
var slope_regrid = slope
    // // Force the next reprojection to aggregate instead of resampling. (unnecesarry here)
    // .reduceResolution({
    //   reducer: ee.Reducer.mean(),
    //   maxPixels: 196530// 1024
    // })
    // // Request the data at the scale and projection of the MODIS image.
    // .reproject({
    //   crs: LSProjection
    // });
var aspect_regrid = aspect
    // // Force the next reprojection to aggregate instead of resampling. (unnecesarry here)
    // .reduceResolution({
    //   reducer: ee.Reducer.mean(),
    //   maxPixels: 196530// 1024
    // })
    // // Request the data at the scale and projection of the MODIS image.
    // .reproject({
    //   crs: LSProjection
    // });
    
///////////////// FUNCTION DEFINITIONS /////////////////////////////////////////////////

/// TOPOGRAPHIC CORRECTION ///

var topo_corr = function(image){
// get solar azimuth and zenith angles from landsat metadata
 var solar_azimuth = ee.Number(LS8_all.first().get('SUN_AZIMUTH')).multiply(Math.PI/180);  
 var solar_zenith = ee.Number(90).subtract((LS8_all.first().get('SUN_ELEVATION'))).multiply(Math.PI/180);
 
 // compute the cos i value for image following Wu et al. 2016 
 var cosiI = ee.Image.constant(solar_zenith).cos().multiply(slope_regrid.cos())
 var cosiII =  ee.Image.constant(solar_zenith).sin().multiply(slope_regrid.sin())
 var cosiIII = (ee.Image.constant(solar_azimuth).subtract(aspect_regrid)).cos()
 var cosiIV = cosiII.multiply(cosiIII)
 var cosi = cosiI.add(cosiIV)
 
 // the Minnaert constant 
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
  
 // print(ui.Chart.array.values(ee.Array(bss), 0, means));
  
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
  
return snow_map.set('threshold',threshold) ; 

}; 


//////////////////////////////////////////////////////////////////////////////////////////////////
// Apply functions to collection
print('ls8',LS8)
// topographical corrections        ////////// !!!!! turn on or off depending on area
var LS_topo_corrected = LS8.map(rad_corr)//.map(topo_corr)
print('After topo corr:',LS_topo_corrected)

// clip the collection to area of interest and band of interest 
var LS_collection= LS_topo_corrected.select('SR_B5') //.map(rad_corr); 

print('ready to run snowmapper:', LS_collection)

// compute the snowmaps 
var snowmaps = LS_collection.map(snowmapper); 

print('Snowmaps:',snowmaps)

print(snowmaps.first())

// also get a list of the tresholds 
var Otsu_thresholds = snowmaps.aggregate_array('threshold') //.reduceColumns(ee.Reducer.toList(), ["snowmap"]).get('list') ;
print('Otsu_thresholds:',Otsu_thresholds);


/////////////////////////////////////////////////////////////////////////////////////////////////
// COMPUTING AAR: 

// Reduce the region. The region parameter is the Feature geometry.
// here we count the number of pixels that are non NAN, we can divide this by the total nb of pixels to get AAR
var AAR = function(image){
    var area = ee.Feature(image.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: GLIMS.first().geometry(),
    scale: 30,
    maxPixels: 1e9}));
    
  return image.set('area', area);
  };
var areas = snowmaps.map(AAR) ; 

print('areas:',areas)
// now et a list of the areas.
var AARs = areas.aggregate_array('area')//.reduceColumns(ee.Reducer.toList(), ["snowmap"]).get('list') ;
print('AARs', AARs);


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

//load external export function 
// var batch = require('users/fitoprincipe/geetools:batch'); 

// batch.Download.ImageCollection.toDrive(snowmaps, 'GIS_final_proj', 
//   {scale: 30, 
//   region:GLIMS, 
//  type: 'float'}); 

// batch.Download.ImageCollection.toDrive(cloud_cover, 'GIS_final_proj_cloudmaps', 
//   {scale: 30, 
//   region:GLIMS, 
//   type: 'float'}); 


// // export the glims polygons 
// var GLIMScollection = GLIMS.reduceToVectors({
//   crs:
// })]) ;

// export the glacier outlines for which the analysis was conducted 
Export.table.toDrive({collection:GLIMS,
  description:'glims_shapefile',
  fileFormat:'kml'
}); 



// visualise the output //////////////////////////////////////////////////////

Map.centerObject(GLIMS,10)

// visaluisation params for glacier outlines
var visGLIMS = {
  palette: ['blue'],
  min: -1.0,
  max: 1.0,
  opacity: 0.4,
};

// visualisation parameters for snow cover 
var visSNOWMAP = {
  palette: ['red'],
  min:-1,
  max:10,
  opacity: 0.4,
};

// add both layers to maps
Map.addLayer(imageGLIMS, visGLIMS, 'GLIMS');
Map.addLayer(snowmaps.first().clip(GLIMS),visSNOWMAP,'snowmap')

//Display B5 clipped to the glaciers: 
var visLS8 = {bands:['SR_B5'],min: -0.1, max: 1}
// add B5 to  to map 
Map.addLayer(LS_collection.first().select(['SR_B5']).clip(GLIMS), visLS8, 'glacier Compsite')

// Make an image histogram 
var chart =
    ui.Chart.image.histogram({image: LS_collection.mean().select("SR_B5"), region: GLIMS, scale: 500})
        .setSeriesNames(['NIR (0.6 to 0.7 µm)'])
        .setOptions({
          title: 'Landsat-8 B5 (NIR: 0.6 to 0.7 µm)',
          hAxis: {
            title: 'Surface Refletance',
            titleTextStyle: {italic: false, bold: true},
          },
          vAxis:
              {title: 'Pixel Count', titleTextStyle: {italic: false, bold: true}},
          colors: ['1d6b99', 'f0af07']
        });
print(chart);

