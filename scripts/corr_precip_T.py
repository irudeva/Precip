from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as datetime  # Python standard library datetime  module
from scipy import interpolate, stats


# pi = np.pi
# dtr=pi/180
# rtd=180/pi
# radius=6.371e6 #radius of sphere having same volume as Earth (m)
# e_omega=7.292e-5 #rotation rate of Earth (rad/s)

#regression test

# X=np.array([0,1,2,3,4,5])
# Y = np.array([np.NaN,4, 5, 10, 2, 5])
# mask = np.isfinite([X, Y]).all(axis=0)
# Xclear = X[mask]
# Yclear = Y[mask]
# print Xclear
# print Yclear
# a, b, r, pval, err = stats.linregress(Xclear, Yclear)
# print r, pval
# print a*Xclear+b

#end regression test




precipin = "../../../DATA/GPCP/precip.mon.mean.nc"
Tin = "../../../DATA/ERAint/air.erain.mon.mean.nc"


yrs=([1979,2016])
ssn   = np.array(['DJF','MAM','JJA','SON'])
for issn in ssn :
    if issn == 'DJF':
        smon = np.array([12, 1, 2])
    elif issn == 'JJA':
        smon = np.array([6,7,8])
    elif issn == 'MAM':
        smon = np.array([3,4,5])
    elif issn == 'SON':
        smon = np.array([9,10,11])

# read precip
precipnam=('lon','lat','time','precip')

print 'precip file in=',precipin

ncP = Dataset(precipin, 'r')

for i, var in enumerate(precipnam):
    print 'precip var: ',var, i
    if ncP.variables[precipnam[i]].name != var:
        print "Variables don't agree", var, nc.variables[varnam[v]].name, v
        exit()

lonp = ncP.variables[precipnam[0]][:]
latp = ncP.variables[precipnam[1]][:]
timep = ncP.variables[precipnam[2]][:]
precip = ncP.variables[precipnam[3]][:]



# read Tin
Tnam=('longitude','latitude','time','t2m')

print 'T file in=',Tin

ncT = Dataset(Tin, 'r')

for i, var in enumerate(Tnam):
    print 't var: ',var, i
    if ncT.variables[Tnam[i]].name != var:
        print "Variables don't agree", var, nc.variables[varnam[v]].name, v
        exit()

lont = ncT.variables[Tnam[0]][:]
latt = ncT.variables[Tnam[1]][:]
timet = ncT.variables[Tnam[2]][:]
if ncT.variables[Tnam[3]].units == "K" :
    t2m = ncT.variables[Tnam[3]][:]-273.15
else:
    t2m = ncT.variables[Tnam[3]][:]

print("Data uploaded")
# end read

#time conversion

dt_timep = [datetime.date(1800, 1, 1) + datetime.timedelta(days=int(t))\
           for t in timep]
dt_timet = [datetime.date(1900, 1, 1) + datetime.timedelta(hours=int(t))\
           for t in timet]



nump = np.zeros([ssn.size, yrs[1]-yrs[0]+1],dtype=int)
numt = np.zeros_like(nump)
precip_ssn = np.zeros([yrs[1]-yrs[0]+1,ssn.size,np.size(latp),np.size(lonp)])
t2m_ssn = np.zeros([yrs[1]-yrs[0]+1,ssn.size,np.size(latt),np.size(lont)])

print 'precip ssn aver shape: ', precip_ssn.shape
print 't2m ssn aver shape: ',t2m_ssn.shape


for yr in range(1979,2017) :
    for js,issn in enumerate(ssn) :
        if issn == 'DJF':
            smon = np.array([12, 1, 2])
        elif issn == 'JJA':
            smon = np.array([6,7,8])
        elif issn == 'MAM':
            smon = np.array([3,4,5])
        elif issn == 'SON':
            smon = np.array([9,10,11])
        # print yr,issn,smon
        for imon in smon :
            iyr = yr
            if imon == 12:
                iyr = yr-1
            for it,t in enumerate(dt_timep) :
                if t == datetime.date(iyr,imon,1):
                    nump[js,yr-yrs[0]]=nump[js,yr-yrs[0]]+1
                    precip_ssn[yr-yrs[0],js,:,:]=precip[it,:,:]+precip_ssn[yr-yrs[0],js,:,:]
            for it,t in enumerate(dt_timet) :
                if t == datetime.date(iyr,imon,1):
                    numt[js,yr-yrs[0]]=numt[js,yr-yrs[0]]+1
                    t2m_ssn[yr-yrs[0],js,:,:]=t2m[it,:,:]+t2m_ssn[yr-yrs[0],js,:,:]

# for yr in range(1979,2017) :
#     for js,issn in enumerate(ssn) :
#         print issn,yr
#         print nump[js,yr-yrs[0]]
#         print numt[js,yr-yrs[0]]
#
#
# print dt_timep
# print dt_timet

for  yr in range(1979,2017) :
    for js,issn in enumerate(ssn) :
        if nump[js,yr-yrs[0]]==3 &   numt[js,yr-yrs[0]]==3 :
            precip_ssn[yr-yrs[0],js,:,:]=precip_ssn[yr-yrs[0],js,:,:]/3
            t2m_ssn[yr-yrs[0],js,:,:]=t2m_ssn[yr-yrs[0],js,:,:]/3
        else :
            print 'deleted season: ',yr, issn
            precip_ssn[yr-yrs[0],js,:,:]=np.nan
            t2m_ssn[yr-yrs[0],js,:,:]=np.nan

##---Interpolation-----------------------------------------------------------------
print "\nInterpolation"

t2m_like_precip = np.zeros_like(precip_ssn)
for  yr in range(1979-yrs[0],2017-yrs[0]) :
    for js,issn in enumerate(ssn) :
        print issn, yr+yrs[0]
        if numt[js,yr]==3 :
            t2mint = interpolate.interp2d(lont, latt, t2m_ssn[yr,js,:,:], kind='cubic')
            for y,lat in enumerate(latp) :
                for x,lon in enumerate(lonp):
                    # print 'lon=',x,lon
                    # print 'lat=',y,lat
                    t2m_like_precip[yr,js,y,x]=t2mint(lon,lat)


print "End Interpolation\n"

print "Control netcdf write:\n"


#control Interpolation
# t2m_ssn to netcdf
fssn = '../output/test/t2m.ssn.nc'
print 'first control netcdf ', fssn
ncout_ssn = Dataset(fssn, 'w', format='NETCDF4')
ncout_ssn.description = "TEST seasonal 2m air temp from %s" % (Tin)

varnam = (Tnam[0],Tnam[1],'ssn','yr','t2m_ssn')
dimnam = (Tnam[0],Tnam[1],'ssn','yr','nchar')

ncout_ssn.createDimension(dimnam[0], lont.size)
ncout_ssn.createDimension(dimnam[1], latt.size)
ncout_ssn.createDimension(dimnam[2], ssn.size)
ncout_ssn.createDimension(dimnam[3], yrs[1]-yrs[0]+1)
ncout_ssn.createDimension(dimnam[4], 3)

for n,nv in enumerate(varnam[:2]) :
    ncout_var = ncout_ssn.createVariable(nv,lont.dtype,dimnam[n])
    for ncattr in ncT.variables[nv].ncattrs():
        ncout_var.setncattr(ncattr, ncT.variables[nv].getncattr(ncattr))

ncout_issn = ncout_ssn.createVariable(varnam[2], 'S1',(dimnam[2],dimnam[4]))
str_out = netCDF4.stringtochar(np.array(ssn, 'S3'))
ncout_yr = ncout_ssn.createVariable(varnam[3], 'i4',dimnam[3])

ncout_var = ncout_ssn.createVariable(varnam[4], 'f',dimnam[3::-1])
ncout_var.long_name = 'seasonally averaged T at 2m'
ncout_var.units = 'C'
# ncout_var.scale_factor = varlist["scale"][iv]
# ncout_var.add_offset   = 0.

ncout_ssn.variables[varnam[0]][:] = lont
ncout_ssn.variables[varnam[1]][:] = latt
ncout_issn[:]                     = str_out
ncout_yr[:]                       = range(yrs[0],yrs[1]+1)
ncout_var[:]                      = t2m_ssn[:,:,:,:]

ncout_ssn.close()

# t2m_like_precip (i.e. t2m interpolated onto precip grid) to netcdf
finter = '../output/test/t2m.inter.nc'
print 'first control netcdf ', finter
ncout_inter = Dataset(finter, 'w', format='NETCDF4')

ncout_inter.description = "TEST seasonal 2m air temp from %s" % (Tin)

varnam = (precipnam[0],precipnam[1],'ssn','yr','t2m_inter')
dimnam = (precipnam[0],precipnam[1],'ssn','yr','nchar')

ncout_inter.createDimension(dimnam[0], lonp.size)
ncout_inter.createDimension(dimnam[1], latp.size)
ncout_inter.createDimension(dimnam[2], ssn.size)
ncout_inter.createDimension(dimnam[3], yrs[1]-yrs[0]+1)
ncout_inter.createDimension(dimnam[4], 3)

for n,nv in enumerate(varnam[:2]) :
    ncout_var = ncout_inter.createVariable(nv,lonp.dtype,dimnam[n])
    for ncattr in ncP.variables[nv].ncattrs():
        ncout_var.setncattr(ncattr, ncP.variables[nv].getncattr(ncattr))

ncout_issn = ncout_inter.createVariable(varnam[2], 'S1',(dimnam[2],dimnam[4]))
ncout_yr   = ncout_inter.createVariable(varnam[3], 'i4',dimnam[3])
ncout_var  = ncout_inter.createVariable(varnam[4], 'f',dimnam[3::-1])
ncout_var.long_name = 'interpolated T at 2m of seasonally averaged values'
ncout_var.units     = 'C'
# ncout_var.scale_factor = varlist["scale"][iv]
# ncout_var.add_offset   = 0.


ncout_inter.variables[varnam[0]][:] = lonp
ncout_inter.variables[varnam[1]][:] = latp
ncout_issn[:]                       = str_out
ncout_yr[:]                         = range(yrs[0],yrs[1]+1)
ncout_var[:]                        = t2m_like_precip

ncout_inter.close()

print "End netcdf control writing\n"
print "End Interpolation\n"

print "------Correlation\n"

for issn in range(ssn.size) :
    for y,lat in enumerate(latp) :
        for x,lon in enumerate(lonp) :
          slope, intercept, r_value, p_value, std_err = stats.linregress(precip_ssn[:,issn,y,x],t2m_like_precip[:,issn,y,x])
          print "precip: "
          print precip_ssn[:,issn,y,x]
          print "t2m   : "
          print t2m_like_precip[:,issn,y,x]
          print ssn[issn], lon,lat, r_value, p_value





ncT.close()
ncP.close()
quit()







print lont[5]
print latt[5]
print t2m_ssn[5,5,3,0]

print t2mint(lont[5],latt[5])



quit()

## plot netcdf fields
## compare the t2m_ssn and inter t2m_ssn

## correlation
## write netcdf of corr coef (regression line)
##plot corr coeff for four seasons
