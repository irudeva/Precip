from netCDF4 import Dataset
import numpy as np
import datetime as datetime  # Python standard library datetime  module


# pi = np.pi
# dtr=pi/180
# rtd=180/pi
# radius=6.371e6 #radius of sphere having same volume as Earth (m)
# e_omega=7.292e-5 #rotation rate of Earth (rad/s)



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

print precipin

nc = Dataset(precipin, 'r')

for i, var in enumerate(precipnam):
    print var, i
    if nc.variables[precipnam[i]].name != var:
        print "Variables don't agree", var, nc.variables[varnam[v]].name, v
        exit()

lonp = nc.variables[precipnam[0]][:]
latp = nc.variables[precipnam[1]][:]
timep = nc.variables[precipnam[2]][:]
precip = nc.variables[precipnam[3]][:]



# read Tin
Tnam=('longitude','latitude','time','t2m')

print Tin

nc = Dataset(Tin, 'r')

for i, var in enumerate(Tnam):
    print var, i
    if nc.variables[Tnam[i]].name != var:
        print "Variables don't agree", var, nc.variables[varnam[v]].name, v
        exit()

lont = nc.variables[Tnam[0]][:]
latt = nc.variables[Tnam[1]][:]
timet = nc.variables[Tnam[2]][:]
t2m = nc.variables[Tnam[3]][:]

print("Data uploaded")
# end read

#time conversion

dt_timep = [datetime.date(1800, 1, 1) + datetime.timedelta(days=int(t))\
           for t in timep]
dt_timet = [datetime.date(1900, 1, 1) + datetime.timedelta(hours=int(t))\
           for t in timet]



nump = np.zeros([ssn.size, yrs[1]-yrs[0]+1],dtype=int)
numt = np.zeros_like(nump)
precip_ssn = np.zeros([np.size(latp),np.size(lonp),ssn.size, yrs[1]-yrs[0]+1])
t2m_ssn = np.empty([np.size(latt),np.size(lont),ssn.size, yrs[1]-yrs[0]+1])

print precip_ssn.shape
print t2m_ssn.shape


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
        print yr,issn,smon
        for imon in smon :
            iyr = yr
            if imon == 12:
                iyr = yr-1
            for it,t in enumerate(dt_timep) :
                if t == datetime.date(iyr,imon,1):
                    nump[js,iyr-yrs[0]]=nump[js,iyr-yrs[0]]+1
                    precip_ssn[:,:,js,iyr-yrs[0]]=precip[it,:,:]+precip_ssn[:,:,js,iyr-yrs[0]]
            for it,t in enumerate(dt_timet) :
                if t == datetime.date(iyr,imon,1):
                    numt[js,iyr-yrs[0]]=numt[js,iyr-yrs[0]]+1
                    t2m_ssn[:,:,js,iyr-yrs[0]]=t2m[it,:,:]+t2m_ssn[:,:,js,iyr-yrs[0]]

for yr in range(1979,2017) :
    for js,issn in enumerate(ssn) :
        print issn,yr
        print nump[js,iyr-yrs[0]]
        print numt[js,iyr-yrs[0]]



quit()

precip_ssn[lon,lat,ssn,year]

t2m_ssn()




quit()


    # nt=np.array([0 for i in range(time.size)])
    # i =0
    # for yr in range(1980,2011) :
    #     for m in bgmon :
    #         yr1 = yr
    #         if m == 12:
    #             yr1 = yr-1
    #         for t in dt_time :
    #             if t == datetime.date(yr1,m,1):
    #                 print 'selected time: ', t
    #                 ind = dt_time.index(t)
    #                 nt[i] = ind
    #                 i += 1
    #
    #
    # u = np.average(uwnd[nt[nt>0],:,:],axis=0)
    # v = np.average(vwnd[nt[nt>0],:,:],axis=0)
    #
    #
    # #----Mercator---------------------------------------------------------------------
    # xm=lons*radius*dtr
    # xm360=360*radius*dtr
    # ym=lats+1  #array declaration
    # #ym[1:-2] = lats[1:-2]
    # ym[1:-2]=radius*np.log((1+np.sin(dtr*lats[1:-2]))/np.cos(dtr*lats[1:-2]));
    # ym[0]=float('inf')
    # ym[-1]=ym[0]
    #
    # # dy = np.gradient(ym)
    # # dx = np.gradient(xm)
    #
    # coslat=np.cos(dtr*lats)
    # #coslat[0]=0   # a very small number is used instead
    # #coslat[-1]=0  # ----"""----
    #
    # # velocity in the Mercator projection
    # um=u/coslat[:,None]
    # vm=v/coslat[:,None]
    #
    #
    # #----BetaM---------------------------------------------------------------------
    # print 'Calculate BetaM'
    #
    # cos2=coslat*coslat
    #
    # # for i in range(0,np.size(um,axis=1)) :
    # #  cosuy_np = np.gradient(um[:,i]*cos2,2*pi*radius/360)
    # #  cosuyy_np = np.gradient(um[:,i]/cos2,2*pi*radius/360)
    # cosuy_np = np.zeros_like(um)
    # cosuyy_np = np.zeros_like(um)
    # #print dy
    # #quit()
    # dy = np.gradient(ym)
    # cosuy_np[:,0] = np.gradient(um[:,0]*cos2,dy)
    # cosuyy_np[:,0] = np.gradient(cosuy_np[:,0]/cos2,dy)
    # for j in range(0,np.size(um,axis=0)) :
    #     cosuy_np[j,:] = cosuy_np[j,0]
    #     cosuyy_np[j,:] = cosuyy_np[j,0]
    #
    # tmp = 2*e_omega *cos2/radius
    # BetaM_np = tmp[:,None]-cosuyy_np
    # BetaM = BetaM_np
    #
    # Ks = np.sqrt(BetaM/um)
    #
    #
    #
    # print("Ks done")
    #
    # #---NetCDF write---------------------------------------------------------------
    # print("Start NetCDF writing")
    #
    # ncout = Dataset(fout, 'w', format='NETCDF4')
    # ncout.description = "Ks %s" % (fout)
    #
    # # Using our previous dimension info, we can create the new time dimension
    # # Even though we know the size, we are going to set the size to unknown
    #
    # ncout.createDimension(dimnam[0], lons.size)
    # ncout.createDimension(dimnam[1], lats.size)
    # #ncout.createDimension(dimnam[2], None)
    #
    # for nv in range(0, 2) :
    #     ncout_var = ncout.createVariable(varnam[nv], nc.variables[varnam[nv]].dtype,dimnam[nv])
    #     for ncattr in nc.variables[varnam[nv]].ncattrs():
    #         ncout_var.setncattr(ncattr, nc.variables[varnam[nv]].getncattr(ncattr))
    # #print(nc.variables['latitude'].ncattrs())
    #
    # ncout.variables[dimnam[0]][:] = lons
    # ncout.variables[dimnam[1]][:] = lats
    # #ncout.variables[dimnam[2]][:] = time
    #
    # ncout_Ks = ncout.createVariable('Ks', 'f',(dimnam[1],dimnam[0]))
    # ncout_Ks.long_name = 'Total stationary wavenumber'
    # #Ks_scale = 1.e-7
    # #Ks_add   = 0.
    # #ncout_Ks.scale_factor = Ks_scale
    # #ncout_Ks.add_offset   = Ks_add
    # #ncout_sf.units        = 'm**2 s**-1'
    #
    # #!!!automatically takes scale and offset into account
    # #!!! no need for: ncout_sf[:] = (sf-sf_add)/sf_scale
    # ncout_Ks[:] = Ks*radius
    #
    #
    # nc.close()
    # ncout.close()
