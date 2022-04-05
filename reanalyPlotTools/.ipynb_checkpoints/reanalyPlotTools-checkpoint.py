import glob
import netCDF4 as nc
import numpy as np
import datetime
from datetime import timedelta, date

def sumSWE(SWE,lenx,leny,resolution):
#     SWE = np.sum(SWE,axis=(1,2))*lenx*leny*resolution
    SWE = np.sum(SWE,axis=(1,2))*(resolution**2)
    return SWE
def sumSCA(SCA,lenx,leny,resolution):
#     SCA = np.sum(SCA,axis=(1,2))*lenx*leny*resolution
    SCA = np.sum(SCA,axis=(1,2))*(resolution**2)
    return SCA
def habitat_Apr(SWE,lenx,leny,resolution,density):
    hab = np.sum(np.where(SWE > 0.33,1,0),axis=(0,1))*lenx*leny*resolution
    return hab
def habitat_May(SWE,lenx,leny,resolution,density):
    hab = np.sum(np.where(SWE > 0.20,1,0),axis=(0,1))*lenx*leny*resolution
    return hab
def date_range(start, end):
    delta = end - start 
    days = [start + timedelta(days=i) for i in range(delta.days + 1)]
    return days

class reanalyPlotTools:
    
    def peakSWE_distribution(dat_path,extension,posteriorIdx,resolution):
        peakSWE = []
        fnames = sorted(glob.glob(dat_path+extension))
        for file in fnames:
            dset = nc.Dataset(file)
            SWE = sumSWE(dset['SWE_Post'][:,posteriorIdx,:,:],
                         len(dset['Longitude'][:]),
                         len(dset['Latitude'][:]),
                         resolution)
            peakSWE.append(np.max(SWE))
        return peakSWE
    
    def P_vs_SWEP(dat_path,extensionSWE,extensionForcing,posteriorIdx,peakSWE,resolution):
        PSUM,SWE_P = [],[]
        fnamesSWE = sorted(glob.glob(dat_path+extensionSWE))
        fnamesForce = sorted(glob.glob(dat_path+extensionForcing))
        for fileSWE,sweVal,fileForce in zip(fnamesSWE,peakSWE,fnamesForce):
            dset = nc.Dataset(fileSWE)
            SWEsave = dset['SWE_Post'][:,posteriorIdx,:,:]
            SWE = sumSWE(dset['SWE_Post'][:,posteriorIdx,:,:],
                         len(dset['Longitude'][:]),
                         len(dset['Latitude'][:]),
                         resolution)            
            end_idx = np.where(SWE == sweVal)[0][0]
            SWE = SWE[0:end_idx]
            start_idx = np.subtract(SWE,sweVal*0.05)
            start_idx = np.where(start_idx < 0,9999,start_idx)
            start_idx = np.where(start_idx == min(start_idx))[0][0]
            SWE = np.mean(SWEsave[start_idx:end_idx,:,:],axis=(1,2))
            dset.close()
            dset = nc.Dataset(fileForce)
            P = np.cumsum(dset['PPT_Post'][0:end_idx,:,:],axis=0)
            PSUM.append(np.mean(P[-1,:,:]-P[start_idx,:,:],axis=(0,1)))
            P = np.mean(P,axis=(1,2))/1000
#             P = np.mean(P[start_idx:],axis=(1,2))/1000
#             P = np.where(P < SWE,np.nan,P)
            P[start_idx:] = np.where(P[start_idx:] < SWE,np.nan,P[start_idx:])
            SWE_P.append(np.mean(SWE/P[start_idx:]))
#             SWE_P.append(np.mean(SWE/P))
            dset.close()
        return PSUM,SWE_P
    
    def extractData(dat_path,extensionSWE,posteriorIdx,resolution,dayApr15,dayMay15,density,PSUM,SWE_P,pctile,years):
        SWEsave,SCAsave,ttSave = np.empty((0,360)),np.empty((0,360)),np.empty((0,360))
        hot_cold,wet_dry,habApr15,habMay15 = [],[],[],[]
        fnamesSWE = sorted(glob.glob(dat_path+extensionSWE))
        for file,psum,swe_p,yr in zip(fnamesSWE,PSUM,SWE_P,years):
            dset = nc.Dataset(file)
            SWE = sumSWE(dset['SWE_Post'][:,posteriorIdx,:,:],
                         len(dset['Longitude'][:]),
                         len(dset['Latitude'][:]),
                         resolution)
            SWEsave = np.vstack((SWEsave,SWE[0:360]))
            SCA = sumSCA(dset['SCA_Post'][:,posteriorIdx,:,:],
                         len(dset['Longitude'][:]),
                         len(dset['Latitude'][:]),
                         resolution)
            SCAsave = np.vstack((SCAsave,SCA[0:360]))
            dtt = date_range(date(yr-1,10,1),date(yr,9,30))
            ttSave = np.vstack((ttSave,dtt[0:360]))
            habApr15.append(habitat_Apr(dset['SWE_Post'][dayApr15,posteriorIdx,:,:],
                         len(dset['Longitude'][:]),
                         len(dset['Latitude'][:]),
                         resolution,density))
            habMay15.append(habitat_May(dset['SWE_Post'][dayMay15,posteriorIdx,:,:],
                         len(dset['Longitude'][:]),
                         len(dset['Latitude'][:]),
                         resolution,density))
            if swe_p > np.percentile(SWE_P,50+pctile):
                hot_cold.append(1)
            elif swe_p < np.percentile(SWE_P,50-pctile):
                hot_cold.append(2)
            else:
                hot_cold.append(0)
            if psum > np.percentile(PSUM,50+pctile):
                wet_dry.append(1)
            elif psum < np.percentile(PSUM,50-pctile):
                wet_dry.append(2)
            else:
                wet_dry.append(0)
        return SWEsave,SCAsave,ttSave,habApr15,habMay15,hot_cold,wet_dry
            
            

    
            
            
            
            