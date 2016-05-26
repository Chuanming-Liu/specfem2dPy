import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy.ma as ma
import matplotlib.animation as animation
import obspy
import os

class StaInfo(object):
    """
    An object contains a station information several methods for station related analysis.
    -----------------------------------------------------------------------------------------------------
    General Parameters:
    stacode       - station name
    network      - network
    chan           - channels for analysis
    x, z             - position for station (in meters)
    -----------------------------------------------------------------------------------------------------
    """
    def __init__(self, stacode=None, network='ME2D', x=None, z=None):

        self.stacode=stacode
        self.network=network
        self.x=x
        self.z=z
        return
    
class StaLst(object):
    """
    An object contains a station list(a list of StaInfo object) information several methods for station list related analysis.
        stations: list of StaInfo
    """
    def __init__(self,stations=None):
        self.stations=[]
        if isinstance(stations, StaInfo):
            stations = [stations]
        if stations:
            self.stations.extend(stations)

    def __add__(self, other):
        """
        Add two StaLst with self += other.
        """
        if isinstance(other, StaInfo):
            other = StaLst([other])
        if not isinstance(other, StaLst):
            raise TypeError
        stations = self.stations + other.stations
        return self.__class__(stations=stations)

    def __len__(self):
        """
        Return the number of Traces in the StaLst object.
        """
        return len(self.stations)

    def __getitem__(self, index):
        """
        __getitem__ method of obspy.Stream objects.
        :return: Trace objects
        """
        if isinstance(index, slice):
            return self.__class__(stations=self.stations.__getitem__(index))
        else:
            return self.stations.__getitem__(index)

    def append(self, station):
        """
        Append a single StaInfo object to the current StaLst object.
        """
        if isinstance(station, StaInfo):
            self.stations.append(station)
        else:
            msg = 'Append only supports a single StaInfo object as an argument.'
            raise TypeError(msg)
        return self

    def ReadStaList(self, stafile):
        """
        Read Sation List from a txt file
        stacode network x z
        """
        f = open(stafile, 'r')
        Sta=[]
        for lines in f.readlines():
            lines=lines.split()
            stacode=lines[0]
            network=lines[1]
            x=float(lines[2])
            z=float(lines[3])
            # if len(lines)==5:
            #     try:
            #         ccflag=int(lines[3])
            #         network=lines[4]
            #     except ValueError:
            #         ccflag=int(lines[4])
            #         network=lines[3]
            # if len(lines)==4:
            #     try:
            #         ccflag=int(lines[3])
            #     except ValueError:
            #         network=lines[3]
            netsta=network+'.'+stacode
            if Sta.__contains__(netsta):
                index=Sta.index(netsta)
                if abs(self[index].lon-lon) >0.01 and abs(self[index].lat-lat) >0.01:
                    raise ValueError('Incompatible Station Location:' + netsta+' in Station List!')
                else:
                    print 'Warning: Repeated Station:' +netsta+' in Station List!'
                    continue
            Sta.append(netsta)
            self.append(StaInfo (stacode=stacode, network=network, x=x, z=z ))
            f.close()
        return
    
    def WriteStaList(self, stafile):
        """
        Read Sation List from a txt file
        stacode network x z
        """
        f = open(stafile, 'w')
        for sta in self.stations:
                tempstr='%s %s %1.5e %1.5e 0.0 0.0  \n' %( sta.stacode, sta.network, sta.x, sta.z )
                f.writelines(tempstr)
        f.close()
        return
    
    def AddSingle(self, x, z, stacode=None, network='MEM2D'):
        if stacode==None:
            stacode=str(int(x))+'S'+str(int(z))
        self.append(StaInfo (stacode=stacode, network=network, x=x, z=z ))
        return
    
    def HomoStaLst(self, xmin, xmax, dx, zmin, zmax, dz, network='MEM2D'):
        Nx = int((xmax-xmin)/dx)+1
        Nz = int((zmax-zmin)/dz)+1
        for i in np.arange(Nx):
            for j in np.arange(Nz):
                x=xmin+dx*i
                z=zmin+dz*j
                stacode=str(int(i))+'S'+str(int(j))
                self.append(StaInfo (stacode=stacode, network=network, x=x, z=z ))
        return
    
    def downsample(self, xmin, dx , zmin, dz):
        nSLst=StaLst()
        for sta in self.stations:
            if (sta.x-xmin)%dx !=0 or  (sta.z-zmin)%dz !=0:
                continue
            nSLst.append(sta)
        return nSLst
    
    def GetInventory(self, outfname=None, source='CUCIEI'):
        """
        Get obspy inventory, used for ASDF dataset
        """
        stations=[]
        total_number_of_channels=1
        site=obspy.core.inventory.util.Site(name='01')
        creation_date=obspy.core.utcdatetime.UTCDateTime(0)
        # channels=[]
        for sta in self.stations:
            channels=[obspy.core.inventory.channel.Channel(code='BSH', location_code='01', latitude=sta.z/100000., longitude=sta.x/100000.,
                        elevation=0.0, depth=0.0)]
            station=obspy.core.inventory.station.Station(code=sta.stacode, latitude=sta.z/100000., longitude=sta.x/100000., elevation=0.0,
                    site=site, channels=channels, total_number_of_channels = total_number_of_channels, creation_date = creation_date)
            stations.append(station)
        network=obspy.core.inventory.network.Network(code=sta.network, stations=stations)
        networks=[network]
        inv=obspy.core.inventory.inventory.Inventory(networks=networks, source=source)
        if outfname!=None:
            inv.write(outfname, format='stationxml')
        return inv
    
        
    def SelectStations(self, x=None, z=None, x0=None, z0=None, dist=None, outflag=True, xmin=-1e10, xmax=1e10, zmin=-1e10, zmax=1e10):
        """
        Select a subset of stations from original station list
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        x, z        - if specified, only append stations with x ==
        -----------------------------------------------------------------------------------------------------
        """
        if x==None and z == None and ( x0==None or z0==None or dist==None):
            raise ValueError("At least one of x or y need to be specified!")
        newSLst=StaLst()
        for sta in self.stations:
            if x0!=None or z0!=None or dist!=None:
                stadist=np. sqrt ( (sta.x-x0)**2 +(sta.z-z0)**2 )
                if stadist < dist and outflag==True:
                    continue
                elif stadist > dist and outflag==False:
                    continue
            if x!=None:
                if sta.x==x and sta.z>zmin and sta.z<zmax:
                    newSLst.append(sta)
            if z!=None:
                if sta.z==z and sta.x>xmin and sta.x<xmax:
                    newSLst.append(sta)
        return newSLst
    
    def GetStation(self, datadir):
        newSLst=StaLst()
        for sta in self.stations:
            txtfname = datadir+'/'+sta.network+'.'+sta.stacode+'.'+'BXY.semd'
            if not os.path.isfile(txtfname):
                continue
            newSLst.append(sta)
        return newSLst
    
    
    

