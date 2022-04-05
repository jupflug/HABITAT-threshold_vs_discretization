import geopandas as gpd
from geopandas import GeoSeries
from shapely.geometry import Polygon, Point
import netCDF4 as nc
import glob
import ulmo
import pandas as pd
import numpy as np

def reanalysisPixels(dset,geo_df_list):
    reanalysis_coords = np.empty((0,2))
    for coordinates in geo_df_list:
        diffy = np.abs(dset['Latitude'][:]-coordinates[0])
        diffx = np.abs(dset['Longitude'][:]-coordinates[1])
        reanalysis_coords = np.vstack((reanalysis_coords,[np.where(diffx == diffx.min())[0][0],
                                                          np.where(diffy == diffy.min())[0][0]])) 
    return reanalysis_coords

def fetch(wsdlurl,sitecode,variablecode,start_date,end_date):
    values_df = None
    try:
#       Pull data into dictionary, indexed by station and date
        site_values = ulmo.cuahsi.wof.get_values(wsdlurl,sitecode,variablecode,
                                                start = start_date,end = end_date)
        values_df = pd.DataFrame.from_dict(site_values['values'])
        values_df['datetime'] = pd.to_datetime(values_df['datetime'], utc=True)
        values_df = values_df.set_index('datetime')
        values_df['value'] = pd.to_numeric(values_df['value']).replace(-9999,np.nan)
        values_df = values_df[values_df['quality_control_level_code'] == '1']
#     Throw exception in the case that no good data exists
    except:
        print("unable to fetch %s" % variablecode)
    return values_df

class snotelPlotTools:
    
    def getBoundingBox(dat_path,extension):
        fnames = sorted(glob.glob(dat_path+extension))
        for file in fnames:
            dset = nc.Dataset(file)
            break
        xmin = dset['Longitude'][:].min()
        xmax = dset['Longitude'][:].max()
        ymin = dset['Latitude'][:].min()
        ymax = dset['Latitude'][:].max()
                
        coords = ((xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin))
        pgon = Polygon(coords)

        df = GeoSeries([pgon])
        return df,dset
    
    def pullSNOTEL(wsdlurl,df,dset,variablecode,start_date,end_date):
        value_dict = {}
        title_dict = {}

        sites = ulmo.cuahsi.wof.get_sites(wsdlurl)
        sites_df = pd.DataFrame.from_dict(sites,orient='index')
        sites_df = sites_df.dropna()
        sites_df['geometry'] = [Point(float(loc['longitude']),
                                      float(loc['latitude'])) for loc in sites_df['location']]
        sites_df = sites_df.drop(columns='location')
        sites_df = sites_df.astype({"elevation_m":float})
        sites_gdf = gpd.GeoDataFrame(sites_df, geometry='geometry')
        sites_gdf.crs = {'init':'epsg:4326'}
    
        poly = df.geometry.unary_union
        sites_gdf_filt = sites_gdf[sites_gdf.geometry.intersects(poly)]
        sites_gdf_filt = sites_gdf_filt.assign(siteStr=sites_gdf_filt.index.str[:])
        
        geo_df_list = [[point.xy[1][0], point.xy[0][0]] for point in sites_gdf_filt.geometry]
        
        reanalysis_coords = reanalysisPixels(dset,geo_df_list)
        
        for i,sitecode in enumerate(sites_gdf_filt.index):
            print('%i of %i sites' % (i+1, len(sites_gdf_filt.index)))
            out = fetch(wsdlurl,sitecode,variablecode,start_date,end_date)
            
            if out is not None:
                value_dict[sitecode] = out['value']
                title_dict[sitecode] = sites_gdf_filt['name'][i]
        
        return sites_gdf_filt,geo_df_list,reanalysis_coords,value_dict,title_dict
    
    def overlappingReanalysis(dat_path,extension,value_dict,posteriorIdx,reanalysis_coords):
        fnames = sorted(glob.glob(dat_path+extension))
        keys = value_dict.keys()
        
        for z,file in enumerate(fnames):
            data = nc.Dataset(file)
            for i,key in enumerate(keys):
                print(file,key)
                data_filt = data['SWE_Post'][:,posteriorIdx,reanalysis_coords[i,0],
                                        reanalysis_coords[i,1]]
                if i == 0:
                    dataHolder = data_filt
                else:
                    dataHolder = np.ma.vstack((dataHolder,data_filt))
            if z == 0:
                overlapping_reanalysis = dataHolder
            else:
                overlapping_reanalysis = np.ma.hstack((overlapping_reanalysis,dataHolder))
        return overlapping_reanalysis
        
        

        
