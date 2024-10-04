from astropy.io import fits
import pysynphot as S
import glob
import numpy as np
import matplotlib.pyplot as plt
#components in ir
#primary mirror
#secondary mirror
#mirror pads

def getthermspec(filt,deltat=0,includepastfilter=False):
    bp1=S.ObsBandpass('wfc3,ir,'+filt)
    mode1=bp1.obsmode
    
    waves=np.arange(1798,20003)
    fluxes1=[]
    names1=[]

    thpt1=[]
    thptname=[]

    for i in mode1._throughput_filenames:
        if i!='clear':
            f=fits.open(i)
            thpt1.append(np.interp(waves,f[1].data.WAVELENGTH,f[1].data.THROUGHPUT))
            thptname.append(i.split('/')[-1].replace('_',''))

    thpt1=np.array(thpt1)[::-1]
    productthpt1=np.cumproduct(thpt1,axis=0)
    #nthrpt=[12,12,11,10,9,8,7,6,5,4,3,2,2,2,2]
    #nthrpt=np.array([12,12,11,10,9,8,7,6,5,4,4,3,2,2,2])-1
    #nthrpt=np.array([12,12,11,10,9,8,7,6,5,4,4,3,3,2,1])-1
    #nthrpt=[11,11,10,9,8,7,6,5,4,3,3,2,2,2,1]
    nthrpt=[11,11,10,9,8,7,6,5,4,3,3,2,2,2,2]
    sources=['primary','pads','secondary','pick off','channel select','fold mirror','ir mirror1','ir mirror2','cold mask','refractive correction plate','refractive correction plate','filter','warm ring','ir window','qe']
    throughputs={'primary':0,'secondary':1,'pick off mirror':2,'channel select':3,'fold mirror':4,'ir mirror1':5,'ir mirror2':6,'cold mask':7,'refractive correction plate':8,'filter':9,'ir window':10,'qe':11,'corr':12}
    nthroughputs={'primary':11,
                  'pads':11,
                  'secondary':10,
                  'pick off':9,
                  'channel select':8,
                  'fold mirror':7,
                  'ir mirror1':6,
                  'ir mirror2':5,
                  'cold mask':4,
                  'refractive correction plate':3,
                  'refractive correction plate':3,
                  'filter':2,
                  'warm ring':2,
                  'ir window':1,
                  'qe':1
                  }
                  
    uncorfluxes1=[]
    compcount=0
    temps=[]
    for i in range(len(mode1.thcompnames)):
        #if compcount==9:
        #    compcount+=1
        #    continue

        if compcount>11 and not includepastfilter:
            break
        if mode1.thcompnames[i]!='clear':
            thfiles=glob.glob('/Users/timothycarleton/synphot_dat/cdbs/comp/wfc3/'+mode1.thcompnames[i]+'*_th.fits')
            f=fits.open(thfiles[-1])
            if 'DEFT' in f[1].header:
                fluxes1.append(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT']+deltat)*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY)*productthpt1[nthrpt[compcount]])
                uncorfluxes1.append(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT'])*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY))
                #plt.plot(waves,fluxes1[-1])
                #plt.plot(waves,uncorfluxes1[-1])
                #plt.plot(waves,productthpt1[nthrpt[compcount]]*1E-7)
                #plt.title(f[0].header['COMPNAME'].replace('_',''))
                #plt.xlim(12000,19000)
                #plt.ylim(0,1E-7)
                #plt.show()
                names1.append(f[0].header['COMPNAME'])
                compcount+=1
                temps.append(f[1].header['DEFT'])


    pixscale=.128
    return waves,np.sum(np.array(fluxes1),axis=0)*mode1.primary_area*pixscale*pixscale


def getthermback(filt,deltat=0,includepastfilter=True):
    bp1=S.ObsBandpass('wfc3,ir,'+filt)
    mode1=bp1.obsmode
    
    waves=np.arange(1798,20003)
    fluxes1=[]
    names1=[]

    thpt1=[]
    thptname=[]

    for i in mode1._throughput_filenames:
        if i!='clear':
            f=fits.open(i)
            thpt1.append(np.interp(waves,f[1].data.WAVELENGTH,f[1].data.THROUGHPUT))
            thptname.append(i.split('/')[-1].replace('_',''))

    thpt1=np.array(thpt1)[::-1]
    productthpt1=np.cumproduct(thpt1,axis=0)
    #nthrpt=[12,12,11,10,9,8,7,6,5,4,3,2,2,2,2]
    #nthrpt=np.array([12,12,11,10,9,8,7,6,5,4,4,3,2,2,2])-1
    #nthrpt=np.array([12,12,11,10,9,8,7,6,5,4,4,3,3,2,1])-1
    #nthrpt=[11,11,10,9,8,7,6,5,4,3,3,2,2,2,1]
    nthrpt=[11,11,10,9,8,7,6,5,4,3,3,2,2,2,2]
    sources=['primary','pads','secondary','pick off','channel select','fold mirror','ir mirror1','ir mirror2','cold mask','refractive correction plate','refractive correction plate','filter','warm ring','ir window','qe']
    throughputs={'primary':0,'secondary':1,'pick off mirror':2,'channel select':3,'fold mirror':4,'ir mirror1':5,'ir mirror2':6,'cold mask':7,'refractive correction plate':8,'filter':9,'ir window':10,'qe':11,'corr':12}
    nthroughputs={'primary':11,
                  'pads':11,
                  'secondary':10,
                  'pick off':9,
                  'channel select':8,
                  'fold mirror':7,
                  'ir mirror1':6,
                  'ir mirror2':5,
                  'cold mask':4,
                  'refractive correction plate':3,
                  'refractive correction plate':3,
                  'filter':2,
                  'warm ring':2,
                  'ir window':1,
                  'qe':1
                  }
                  
    uncorfluxes1=[]
    compcount=0
    temps=[]
    
    for i in range(len(mode1.thcompnames)):
        #if compcount==9:
        #    compcount+=1
        #    continue

        if compcount>11 and not includepastfilter:
            break
        if mode1.thcompnames[i]!='clear':
            thfiles=glob.glob('/Users/timothycarleton/synphot_dat/cdbs/comp/wfc3/'+mode1.thcompnames[i]+'*_th.fits')
            f=fits.open(thfiles[-1])
            if 'DEFT' in f[1].header:
                fluxes1.append(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT']+deltat)*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY)*productthpt1[nthrpt[compcount]])
                uncorfluxes1.append(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT'])*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY))
                #plt.plot(waves,fluxes1[-1])
                #plt.plot(waves,uncorfluxes1[-1])
                #plt.plot(waves,productthpt1[nthrpt[compcount]]*1E-7)
                #plt.title(f[0].header['COMPNAME'].replace('_',''))
                #plt.xlim(12000,19000)
                #plt.ylim(0,1E-7)
                #plt.show()
                print(f[0].header['COMPNAME'],f[1].header['DEFT']-273.15)
                names1.append(f[0].header['COMPNAME'])
                compcount+=1
                temps.append(f[1].header['DEFT'])

    print(names1)
    pixscale=.128
    return np.trapz(np.sum(np.array(fluxes1),axis=0),x=waves)*mode1.primary_area*pixscale*pixscale

def getthermback_element(filt,element,deltat=0,includepastfilter=True):
    bp1=S.ObsBandpass('wfc3,ir,'+filt)
    mode1=bp1.obsmode
    
    waves=np.arange(1798,20003)
    fluxes1=[]
    names1=[]

    thpt1=[]
    thptname=[]

    for i in mode1._throughput_filenames:
        if i!='clear':
            f=fits.open(i)
            thpt1.append(np.interp(waves,f[1].data.WAVELENGTH,f[1].data.THROUGHPUT))
            thptname.append(i.split('/')[-1].replace('_',''))

    thpt1=np.array(thpt1)[::-1]
    

    productthpt1=np.cumproduct(thpt1,axis=0)

    
    #nthrpt=[12,12,11,10,9,8,7,6,5,4,3,2,2,2,2]
    #nthrpt=np.array([12,12,11,10,9,8,7,6,5,4,4,3,2,2,2])-1
    #nthrpt=np.array([12,12,11,10,9,8,7,6,5,4,4,3,3,2,1])-1
    #nthrpt=[11,11,10,9,8,7,6,5,4,3,3,2,2,2,1]
    nthrpt=[11,11,10,9,8,7,6,5,4,3,3,2,2,2,2]
    sources=['primary','pads','secondary','pick off','channel select','fold mirror','ir mirror1','ir mirror2','cold mask','refractive correction plate','refractive correction plate','filter','warm ring','ir window','qe']
    throughputs={'primary':0,'secondary':1,'pick off mirror':2,'channel select':3,'fold mirror':4,'ir mirror1':5,'ir mirror2':6,'cold mask':7,'refractive correction plate':8,'filter':9,'ir window':10,'qe':11,'corr':12}
    nthroughputs={'primary':11,
                  'pads':11,
                  'secondary':10,
                  'pick off':9,
                  'channel select':8,
                  'fold mirror':7,
                  'ir mirror1':6,
                  'ir mirror2':5,
                  'cold mask':4,
                  'refractive correction plate':3,
                  'refractive correction plate':3,
                  'filter':2,
                  'warm ring':2,
                  'ir window':1,
                  'qe':1
                  }
                  
    uncorfluxes1=[]
    compcount=0
    temps=[]
    for i in range(len(mode1.thcompnames)):
        #if compcount==9:
        #    compcount+=1
        #    continue

        if compcount>11 and not includepastfilter:
            break
        if mode1.thcompnames[i]!='clear':
            thfiles=glob.glob('/Users/timothycarleton/synphot_dat/cdbs/comp/wfc3/'+mode1.thcompnames[i]+'*_th.fits')
            f=fits.open(thfiles[-1])
            if 'DEFT' in f[1].header:
                print(mode1.thcompnames[i],f[1].header['DEFT'])
                if mode1.thcompnames[i]==element:
                    fluxes1.append(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT']+deltat)*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY)*productthpt1[nthrpt[compcount]])
                    uncorfluxes1.append(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT'])*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY))
                else:
                    fluxes1.append(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT'])*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY)*productthpt1[nthrpt[compcount]])
                    uncorfluxes1.append(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT'])*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY))
                    #plt.plot(waves,fluxes1[-1])
                    #plt.plot(waves,uncorfluxes1[-1])
                    #plt.plot(waves,productthpt1[nthrpt[compcount]]*1E-7)
                    #plt.title(f[0].header['COMPNAME'].replace('_',''))
                    #plt.xlim(12000,19000)
                    #plt.ylim(0,1E-7)
                    #plt.show()
                names1.append(f[0].header['COMPNAME'])
                compcount+=1
                temps.append(f[1].header['DEFT'])
                


    pixscale=.128
    for i in range(len(names1)):
        print(names1[i],np.trapz(fluxes1[i],x=waves)*mode1.primary_area*pixscale*pixscale)
    return np.trapz(np.sum(np.array(fluxes1),axis=0),x=waves)*mode1.primary_area*pixscale*pixscale

def getthermback_justelement(filt,element,deltat=0,includepastfilter=True):
    bp1=S.ObsBandpass('wfc3,ir,'+filt)
    mode1=bp1.obsmode
    
    waves=np.arange(1798,20003)
    fluxes1=[]
    names1=[]

    thpt1=[]
    thptname=[]

    for i in mode1._throughput_filenames:
        if i!='clear':
            f=fits.open(i)
            thpt1.append(np.interp(waves,f[1].data.WAVELENGTH,f[1].data.THROUGHPUT))
            thptname.append(i.split('/')[-1].replace('_',''))

    thpt1=np.array(thpt1)[::-1]
    productthpt1=np.cumproduct(thpt1,axis=0)
    
    #nthrpt=[12,12,11,10,9,8,7,6,5,4,3,2,2,2,2]
    #nthrpt=np.array([12,12,11,10,9,8,7,6,5,4,4,3,2,2,2])-1
    #nthrpt=np.array([12,12,11,10,9,8,7,6,5,4,4,3,3,2,1])-1
    #nthrpt=[11,11,10,9,8,7,6,5,4,3,3,2,2,2,1]
    nthrpt=[11,11,10,9,8,7,6,5,4,3,3,2,2,2,2]
    sources=['primary','pads','secondary','pick off','channel select','fold mirror','ir mirror1','ir mirror2','cold mask','refractive correction plate','refractive correction plate','filter','warm ring','ir window','qe']
    throughputs={'primary':0,'secondary':1,'pick off mirror':2,'channel select':3,'fold mirror':4,'ir mirror1':5,'ir mirror2':6,'cold mask':7,'refractive correction plate':8,'filter':9,'ir window':10,'qe':11,'corr':12}
    nthroughputs={'primary':11,
                  'pads':11,
                  'secondary':10,
                  'pick off':9,
                  'channel select':8,
                  'fold mirror':7,
                  'ir mirror1':6,
                  'ir mirror2':5,
                  'cold mask':4,
                  'refractive correction plate':3,
                  'refractive correction plate':3,
                  'filter':2,
                  'warm ring':2,
                  'ir window':1,
                  'qe':1
                  }
                  
    uncorfluxes1=[]
    compcount=0
    temps=[]
    for i in range(len(mode1.thcompnames)):
        #if compcount==9:
        #    compcount+=1
        #    continue

        if compcount>11 and not includepastfilter:
            break
        if mode1.thcompnames[i]!='clear':
            thfiles=glob.glob('/Users/timothycarleton/synphot_dat/cdbs/comp/wfc3/'+mode1.thcompnames[i]+'*_th.fits')
            f=fits.open(thfiles[-1])
            if 'DEFT' in f[1].header:
                #print(mode1.thcompnames[i],f[1].header['DEFT'])
                if mode1.thcompnames[i]==element:
                    fluxes1.append(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT']+deltat)*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY)*productthpt1[nthrpt[compcount]])
                    uncorfluxes1.append(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT'])*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY))

                    #plt.plot(waves,fluxes1[-1])
                    #plt.plot(waves,uncorfluxes1[-1])
                    #plt.plot(waves,productthpt1[nthrpt[compcount]]*1E-7)
                    #plt.title(f[0].header['COMPNAME'].replace('_',''))
                    #plt.xlim(12000,19000)
                    #plt.ylim(0,1E-7)
                    #plt.show()
                    names1.append(f[0].header['COMPNAME'])
                    temps.append(f[1].header['DEFT'])
                compcount+=1
                


    pixscale=.128
    #for i in range(len(names1)):
        #print(names1[i],np.trapz(fluxes1[i],x=waves)*mode1.primary_area*pixscale*pixscale)
    
    return np.trapz(np.sum(np.array(fluxes1),axis=0),x=waves)*mode1.primary_area*pixscale*pixscale


def getthermspec_element(filt,element,deltat=0,includepastfilter=True):
    #element chose from: ['wfc3_ir_primary', 'wfc3_ir_pads', 'wfc3_ir_secondary', 'wfc3_pom', 'wfc3_ir_csm', 'wfc3_ir_fold', 'wfc3_ir_mir1', 'wfc3_ir_mir2', 'wfc3_ir_mask', 'wfc3_ir_rcp', 'wfc3_ir_rcp', 'wfc3_ir_f160w', 'wfc3_ir_wmring', 'wfc3_ir_win', 'wfc3_ir_qe']
    bp1=S.ObsBandpass('wfc3,ir,'+filt)
    mode1=bp1.obsmode
    
    waves=np.arange(1798,20003)
    fluxes1=[]
    names1=[]

    thpt1=[]
    thptname=[]

    for i in mode1._throughput_filenames:
        if i!='clear':
            f=fits.open(i)
            thpt1.append(np.interp(waves,f[1].data.WAVELENGTH,f[1].data.THROUGHPUT))
            thptname.append(i.split('/')[-1].replace('_',''))

    thpt1=np.array(thpt1)[::-1]
    productthpt1=np.cumproduct(thpt1,axis=0)
    #nthrpt=[12,12,11,10,9,8,7,6,5,4,3,2,2,2,2]
    #nthrpt=np.array([12,12,11,10,9,8,7,6,5,4,4,3,2,2,2])-1
    #nthrpt=np.array([12,12,11,10,9,8,7,6,5,4,4,3,3,2,1])-1
    #nthrpt=[11,11,10,9,8,7,6,5,4,3,3,2,2,2,1]
    nthrpt=[11,11,10,9,8,7,6,5,4,3,3,2,2,2,2]
    sources=['primary','pads','secondary','pick off','channel select','fold mirror','ir mirror1','ir mirror2','cold mask','refractive correction plate','refractive correction plate','filter','warm ring','ir window','qe']
    throughputs={'primary':0,'secondary':1,'pick off mirror':2,'channel select':3,'fold mirror':4,'ir mirror1':5,'ir mirror2':6,'cold mask':7,'refractive correction plate':8,'filter':9,'ir window':10,'qe':11,'corr':12}
    nthroughputs={'primary':11,
                  'pads':11,
                  'secondary':10,
                  'pick off':9,
                  'channel select':8,
                  'fold mirror':7,
                  'ir mirror1':6,
                  'ir mirror2':5,
                  'cold mask':4,
                  'refractive correction plate':3,
                  'refractive correction plate':3,
                  'filter':2,
                  'warm ring':2,
                  'ir window':1,
                  'qe':1
                  }
                  
    uncorfluxes1=[]
    compcount=0
    temps=[]

    fluxes=np.zeros_like(waves).astype(np.float)
    
    for i in range(len(mode1.thcompnames)):
        #if compcount==9:
        #    compcount+=1
        #    continue

        if compcount>11 and not includepastfilter:
            break
        if mode1.thcompnames[i]!='clear':
            thfiles=glob.glob('/Users/timothycarleton/synphot_dat/cdbs/comp/wfc3/'+mode1.thcompnames[i]+'*_th.fits')
            f=fits.open(thfiles[-1])
            if 'DEFT' in f[1].header:
                fluxes1.append(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT']+deltat)*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY)*productthpt1[nthrpt[compcount]])
                uncorfluxes1.append(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT'])*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY))
                fluxes+=np.array(S.planck.bb_photlam_arcsec(waves,f[1].header['DEFT']+deltat)*f[1].header['BEAMFILL']*np.interp(waves,f[1].data.WAVELENGTH,f[1].data.EMISSIVITY)*productthpt1[nthrpt[compcount]]).astype(np.float)
                #plt.plot(waves,fluxes1[-1])
                #plt.plot(waves,uncorfluxes1[-1])
                #plt.plot(waves,productthpt1[nthrpt[compcount]]*1E-7)
                #plt.title(f[0].header['COMPNAME'].replace('_',''))
                #plt.xlim(12000,19000)
                #plt.ylim(0,1E-7)
                #plt.show()
                names1.append(f[0].header['COMPNAME'])
                compcount+=1
                temps.append(f[1].header['DEFT'])


    w=np.where(np.array(names1)==element)[0]
    pixscale=.128
    return waves,np.array(fluxes1)[w]*mode1.primary_area*pixscale*pixscale
