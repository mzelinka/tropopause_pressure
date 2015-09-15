###########################################################################
def get_tropopause_pressure(mo,exp,rip,fprov,realm,tres,tab):

    # Perform tropopause pressure calculation
    # generated via: f2py -c -m tropo tropo.f90 
    # tropo.f90 came from http://www.inscc.utah.edu/~reichler/research/projects/TROPO/code.txt

    import tropo  

    f=PMC_utils.fopen2(mo,exp,rip,'ta',fprov,realm,tres,tab)
    tslc1=slice(0,150*12) #default=include 1st pt, exclude last.
    v_fut=f('ta',time=tslc1,level=(botp,topp,'ccn')) #dims=[time,lev?,lat,lon]
    Ta0=v_fut.regrid(kern_grid,regridTool="esmf")     
    f.close()
    del v_fut

    PP=kern_plev/100.
    plimu=45000
    pliml=7500
    plimlex=7500
    dofill=0
    # fortran wants this in lon,lat,ctp but lets swap time with lon
    tropp=MV.zeros((Ta0[:,0,:].shape))
    tropp=MV.masked_where(tropp==0,tropp)
    for LO in range(len(kern_lons)):
        temp=Ta0[:,:,:,LO](order=[0,2,1]) 
        tp,tperr=tropo.tropo(temp,PP,plimu,pliml,plimlex,dofill)

        tropp[:,:,LO]=tp/100. # convert to hPa

    tropp=MV.masked_where(tropp>plimu/100.,tropp)
    tropp=MV.masked_where(tropp<0.,tropp)

    return tropp
