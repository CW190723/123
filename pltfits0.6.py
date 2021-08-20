#!usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import re
import matplotlib.ticker as ticker
import pandas as pd
import math
import csv
import os
from astropy.wcs import WCS
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from scipy import signal

#----------------------------------------------------------------choisemod-----------------------------------------------------------------:
choisemod = 0.2                        #0,0.1,0.2,1,2,3,5,6


#------------------------------------------------------------------------------------------------------------------------------------------:
#------------------------------------------------------------------common------------------------------------------------------------------:
if choisemod == 1 or choisemod == 2 or choisemod == 3:
    fitspath = '/home/cuiqifan/Desktop/jiaocha/fits/'
    filename = 'hlsp_frontier_hst_acs-60mas_abell2744_f814w_v1.0-epoch2_drz.fits'
    n = 0                            #hdu[n].data
    pixscale = 0.06
    vmin = 18
    vmax = 28
    title = 'a'
    xname = 'Dc(pixel)'
    yname = 'Dc(pixel)'

#------------------------------------------------------------------model1------------------------------------------------------------------:
if choisemod == 1:
    Taxesset=False
    xmin=1000
    xmax=2000
    ymin=1000
    ymax=2000

#------------------------------------------------------------------model2------------------------------------------------------------------:
if choisemod == 2:
    csvpath1 = '/home/cuiqifan/Desktop/jiaocha/catalog/'
    csvname1 = 'hlsp_frontier_hst_acs-60mas_abell2744-hffpar_f814w_v1.0_comb.cat'
#    csvpath2 = 
#    csvname2 = 
    xytype1 = True                   #if true, it will use x and y to plot ellipse; if false, it will use ra and dec to plot ellipse.
    selectname1 = ['NUMBER','X_IMAGE','Y_IMAGE','A_IMAGE','B_IMAGE','KRON_RADIUS','THETA_IMAGE']
#    xythpe2 = 
#    selectname2 = []
    Taxesset=False
    xmin=1000
    xmax=2000
    ymin=1000
    ymax=2000

#------------------------------------------------------------------model3------------------------------------------------------------------:
if choisemod == 3:
    csvpath = '/home/cuiqifan/Desktop/jiaocha/catalog/'
    csvname = 'hlsp_frontier_hst_acs-60mas_abell2744-hffpar_f814w_v1.0_comb.cat'
    xytype = True
    selectname = ['NUMBER','X_IMAGE','Y_IMAGE','A_IMAGE','B_IMAGE','KRON_RADIUS','THETA_IMAGE']
    numrows = 2                      #numrows * numcols is pictrue's number of one figure
    numcols = 3
    kpcpicscale = False              #if true, it will use (picscale)kpc to cut picture.(if you use it, you must give 'z_peak' in selectname.)
    picscale = 25                    #for exmple, if you want to have a 50*50 picture scale, you can write 25

#------------------------------------------------------------------model5------------------------------------------------------------------:
#save ds9 parameter
if choisemod == 5:
    csvpath = ''
    csvname = ''
    savetype = 'ellipse'             #ellipse,box,circle
    selectname = ['OBJ_ID','X','Y','SEMI_MAJOR','SEMI_MINOR','POSITION_ANGLE']
#ellipse is ['ID','X','Y','changzhou','duanzhou','kron','POSITION_ANGLE']; box is ['ID','X','Y','xmin','xmax','ymin','ymax']
#circle is ['ID','X','Y','radius']

#------------------------------------------------------------------model6------------------------------------------------------------------:
if choisemod == 6:
    fitspath1 = '/home/cuiqifan/Documents/kerneltest/'
    fitsname1 = 'gs_psf_f814w.fits'  #be convolve file
    fitspath2 = '/home/cuiqifan/Documents/kerneltest/'
    fitsname2 = 'gs_psf_f814w_f160w.fits'                                #convolve kernel
    n = 0                                                                #hdu
    kn = 0                                                               #kernel hdu

#----------------------------------------------------------------kpctopix(0)---------------------------------------------------------------:
if choisemod == 0:
    z_peak = 0.308
    pixscale = 0.03

#---------------------------------------------------------------savecsv(0.1)---------------------------------------------------------------:
if choisemod == 0.1:
    name = ''
    colname = []
    data = 1

#---------------------------------------------------------------savecsv(0.2)---------------------------------------------------------------:
if choisemod == 0.2:
    fitspath = '/home/cuiqifan/Desktop/jiaocha/fits/'
    fitsname = 'hlsp_frontier_hst_acs-60mas_abell2744-hffpar_f814w_v1.0-epoch1_drz.fits'
    hdu = 0    #0,1,2et
    mode = 'img' #'img','wcs'
    center = '2000,2000' #'XXXX,XXXX'
    width = 200

#---------------------------------------------------------------savecsv(0.3)---------------------------------------------------------------:
if choisemod == 0.3:
    data=[[[hdr0],[img0]],[[hdr1],[img1]],[[hdr2],[img2]]]

#-----------------------------------------------------------MJy/srtoJy/pixel(0.4)----------------------------------------------------------:
if choisemod == 0.4:
    pixscal = 0.6
    SB = 1.0

#------------------------------------------------------------------------------------------------------------------------------------------:
#----------------------------------------------------------------selectdata----------------------------------------------------------------:
openselectdata=False
#def selectdata(data_new):
#    for i in data_new:
#        if i[1]>0:
#            data_new=data_new


#------------------------------------------------------------------------------------------------------------------------------------------:
#-------------------------------------------------------------------code-------------------------------------------------------------------:
if choisemod == 1 or choisemod == 2 or choisemod == 3:
    hdu = fits.open(fitspath+filename)
    img=hdu[n].data
    try:
        zp=re.findall(r'f\d\d\d',filename)
    except:
        print('you need give zp a filter value such as f160 or write it into your filename of fits.')
    zpall={'f105':'26.2687','f125':'26.2303','f140':'26.4524','f160':'25.9463','f435':'25.665','f606':'26.493','f814':'25.947'}
    img[img<0]=1e-7
    img[img==0]=1e-10
    img = -2.5 * (np.log10(img) - np.log10(pixscale**2)) + float(zpall[zp[0]])
    
    Dim1 = hdu[n].header['NAXIS1']
    Dim2 = hdu[n].header['NAXIS2']

#---------------------------------------------------------------kpctopix(0)---------------------------------------------------------------:
def kpctopix(z_peak,pixscale):
    cosmo = FlatLambdaCDM(H0=70,Om0=0.3)
    DC = cosmo.comoving_distance(z_peak).value
    LD = cosmo.luminosity_distance(z_peak).value
    DA = DC/(1+z_peak)
    arcsec_per1kpc = 1/(DA*math.pi*1000.0/(3600.0*180.0))
    global pix_per1kpc
    pix_per1kpc = arcsec_per1kpc/pixscale

#--------------------------------------------------------------savecsv(0.1)---------------------------------------------------------------:
def savecsv(name,colname,data):
    csvFile = open(name,'w')
    writer = csv.writer(csvFile)
    writer.writerow(colname)
    writer.writerows(data)
    csvFile.close()

#--------------------------------------------------------------cutfits(0.2)---------------------------------------------------------------:
if choisemod == 0.2:
    def cutfits(hdu,mode,center,width):
        os.system('astcrop --hdu={} --mode={} --center={} --width={} {}'.format(hdu,mode,center,width,fitspath+fitsname))

#-------------------------------------------------------------savefits(0.3)---------------------------------------------------------------:
def savefits(*data):
    if len(data)>0:
        hdr0=data[0][0]
        img0=data[0][1]
        grey0=fits.PrimaryHDU(img0,hdr0)
    if len(data)>1:
        hdr1=data[1][0]
        img1=data[1][1]
        grey1=fits.ImageHDU(img1,hdr1)
    if len(data)>2:
        hdr2=data[2][0]
        img2=data[2][1]
        grey2=fits.ImageHDU(img2,hdr2)
    if len(data)>3:
        hdr3=data[3][0]
        img3=data[3][1]
        grey3=fits.ImageHDU(img3,hdr3)
    
    if len(data)==1:
        greyHDU=fits.HDUList([grey0])
    if len(data)==2:
        greyHDU=fits.HDUList([grey0,grey1])
    if len(data)==3:
        greyHDU=fits.HDUList([grey0,grey1,grey2])
    if len(data)==4:
        greyHDU=fits.HDUList([grey0,grey1,grey2,grey3])
    greyHDU.writeto('newfits2.fits')

#---------------------------------------------------------MJy/srtoJy/pixel(0.4)-----------------------------------------------------------:
if choisemod == 0.4:
    def MJysr_to_Jypix(pixscal,SB):
        C = 4 * math.pi * (pixscal/3600) ** 2 / 41252
        F = SB * C
        return F

#----------------------------------------------------------------pltfits1-----------------------------------------------------------------:
if choisemod == 1:
    def pltfits1(Taxesset=False):
        fig,axes=plt.subplots(ncols=1, nrows=1, figsize=(30,20*Dim1/Dim2))
        plt.subplots_adjust(top=0.9,bottom=0.1,left=0.07,right=0.95, wspace=0.05)
        axes.set_title('{}'.format(title), fontsize=60)
        axes.set_xlabel('{}'.format(xname), fontsize=50)
        axes.set_ylabel('{}'.format(yname), fontsize=50)
        axes.xaxis.set_major_locator(ticker.NullLocator())
        axes.yaxis.set_major_locator(ticker.NullLocator())
        axes.imshow(img, cmap='gray', origin='lower',aspect='equal',vmin=vmin,vmax=vmax)
        if Taxesset==True:
            axes.set_xlim([xmin,xmax])
            axes.set_ylim([ymin,ymax])
        print(hdu.info())
        plt.savefig('{}.png'.format(filename[:-5]))

#----------------------------------------------------------------pltfits2-----------------------------------------------------------------:
if choisemod == 2:
    def pltfits2(csvname1,csvname2=False,Taxesset=False):
        data = pd.read_csv(csvpath1+csvname1)
        data_new = np.array(data.loc[:,selectname1])
        if openselectdata == True:
            selectdata(data_new)
        if csvname2 != False:
            data2 = pd.read_csv(csvpath2+csvname2)
            data_new2 = np.array(data2.loc[:,selectname2])
            if openselectdata == True:
                selectdata(data_new2)
        
        wcs=WCS(filename)
        fig,axes=plt.subplots(ncols=1, nrows=1, figsize=(30,20*Dim1/Dim2))
        plt.subplots_adjust(top=0.92,bottom=0.07,left=0.07,right=0.95, wspace=0.05)
        hb1=axes.imshow(img, cmap='gray', origin='lower',aspect='equal',vmin=vmin,vmax=vmax)
        for l in data_new:
            if xytype1 == False:
                x,y = wcs.all_world2pix(j[1],j[2],0)
            if xytype1 == True:
                x,y = l[1],l[2]
            e1 = Ellipse(xy = (x-1,y-1), width = 2*l[3]*l[5], height = 2*l[4]*l[5], angle=l[6],fill=False,color='r')
            axes.add_artist(e1)
        if csvname2 != False:
            for l in data_new2:
                if xytype2 == False:
                    x,y = wcs.all_world2pix(j[1],j[2],0)
                if xytype2 == True:
                    x,y = l[1],l[2]
                e2 = Ellipse(xy = (x-1,y-1), width = 2*l[3]*l[5], height = 2*l[4]*l[5], angle=l[6],fill=False,color='y')
                axes.add_artist(e2)
        axes.set_title('{}'.format(title), fontsize=60)
        axes.set_xlabel('{}'.format(xname), fontsize=40)
        axes.set_ylabel('{}'.format(yname), fontsize=40)
        axes.tick_params(labelsize=30)
        axes.imshow(img, cmap='gray', origin='lower',aspect='equal',vmin=vmin,vmax=vmax)
        if Taxesset==True:
            axes.set_xlim([xmin,xmax])
            axes.set_ylim([ymin,ymax])
            plt.tick_params(labelsize=30)
        cb = plt.colorbar(hb1)
        cb.set_label('SB{}'.format(zp),size=40)
        cb.ax.tick_params(labelsize=30)
        plt.savefig('{}.png'.format(filename[:-5]))

#----------------------------------------------------------------pltfits3-----------------------------------------------------------------:
if choisemod == 3:
    def pltfits3(csvname,xytype,numrows,numcols,picscale):
        number = numrows * numcols
        data = pd.read_csv(csvpath+csvname)
        data_new = np.array(data.loc[:,selectname])
        wcs=WCS(filename)
        if openselectdata == True:
            selectdata(data_new)
        with PdfPages('{}.pdf'.format(filename[:-5])) as pdf:
            for i in range(math.ceil(len(data_new)/number)):
                fig,axes=plt.subplots(nrows=numrows,ncols=numcols,figsize=(10*numcols,10*numrows))
                plt.subplots_adjust(top=0.92,bottom=0.07,left=0.07,right=0.95, wspace=0.2)
                fig.suptitle('{}'.format(title),fontsize=40)
                t = 0
                for j in data_new[i*number:(i+1)*number]:
                    if xytype == False:
                        x,y = wcs.all_world2pix(j[1],j[2],0)
                    if xytype == True:
                        x,y = j[1],j[2]
                    if kpcpicscale == True:
                        try:
                            kpctopix(i[7],pixscale)
                            picscale = int(pix_per1kpc/2)
                        except:
                            print('you must give the z_peak!')
                    a=img[int(y)-picscale-1:int(y)+picscale,int(x)-picscale-1:int(x)+picscale]
                    axes.ravel()[t].imshow(a, cmap='gray', origin='lower',aspect='equal',vmin=vmin,vmax=vmax)
                    e = Ellipse(xy = (picscale+(x-int(x)),picscale+(y-int(y))), width = 2*j[3]*j[5], height = 2*j[4]*j[5], angle=j[6],fill=False,color='r')
                    axes.ravel()[t].add_artist(e)
                    axes.ravel()[t].set_title('id={:.0f},pos={:.2f},{:.2f}'.format(j[0],j[1],j[2]),fontsize=30)
                    axes.ravel()[t].set_xlabel('{}'.format(xname),fontsize=30)
                    axes.ravel()[t].set_ylabel('{}'.format(yname),fontsize=30)
                    axes.ravel()[t].tick_params(labelsize=20)
                    t = t+1
                pdf.savefig()
                plt.close()

#----------------------------------------------------------------pltfits4-----------------------------------------------------------------:
if choisemod == 4:
    def searchfilename(wht=False):
        name = ['abell2744','abells1063','abell370','macs0416','macs0717','macs1149','hudfp2','hudfp3','hudfp4','xdf']
        epoch160 = [2,1,1,1,1,2]
        epoch814 = [1,2,2,2,2,1]
        mas = 30
        global hudfp,xdf,hff
        hudfp,xdf,hff = [],[],[]
        for i in range(len(name)):
            if 'hudfp' in name[i]:
                path = '/home/cuiqifan/data/HUDF/{}/hlsp_hlf_hst_{}-{}mas_{}_{}'.format(name[i].upper(),{},mas,name[i],{})
                filename160 = path.format('wfc3','f160w_v2.0_sci.fits')
                filename814 = path.format('acs','f814w_v2.0_sci.fits')
                if wht == True:
                    filenamewht814 = path.format('acs','f814w_v2.0_wht.fits')
                    hudfp.append([filename160,filename814,filenamewht814])
                elif wht == False:
                    hudfp.append([filename160,filename814])
            elif 'xdf' in name[i]:
                path = '/home/cuiqifan/data/HUDF/{}/hlsp_xdf_hst_{}-{}mas_hudf_{}'.format(name[i].upper(),mas,{},{})
                filename160 = path.format('wfc3ir','f160w_v1_sci.fits')
                filename814 = path.format('acswfc','f814w_v1_sci.fits')
                if wht == True:
                    filenamewht814 = path.format('acswfc','f814w_v1_wht.fits')
                    xdf.append([filename160,filename814,filenamewht814])
                elif wht == False:
                    xdf.append([filename160,filename814])
            else:
                path = '/home/cuiqifan/data/{}/epoch{}/{}/hlsp_frontier_hst_{}-{}mas_{}-hffpar_{}_v1.0-epoch{}_{}'.format(name[i],{},{},{},mas,name[i],{},{},{})
                filename160 = path.format(epoch160[i],'f160w','wfc3','f160w',epoch160[i],'drz.fits')
                filename814 = path.format(epoch814[i],'f814w','acs','f814w',epoch814[i],'drz.fits')
                if wht == True:
                    filenamewht814 = path.format(epoch814[i],'f814w','acs','f814w',epoch814[i],'wht.fits')
                    hff.append([filename160,filename814,filenamewht814])
                elif wht == False:
                    hff.append([filename160,filename814])
        
#-----------------------------------------------------------------model5------------------------------------------------------------------:
#save ds9 parameter
if choisemod == 5:
    def pltfits5(csvpath,csvname,selectname,savetype):
        data = pd.read_csv(csvpath+csvname)
        data_new = np.array(data.loc[:,selectname])
        
        str2 = []
        if savetype == 'ellipse':
            for i in data_new:
                str1 = 'ellipse('+str(i[1])+' '+str(i[2])+' '+str(i[3]*i[5])+' '+str(i[4]*i[5])+' '+str(i[6])+')'
                str2.append(str1)
        elif savetype == 'box':
            for i in data_new:
                str1 = 'box('+str(i[1])+' '+str(i[2])+' '+str(i[4]-i[3])+' '+str(i[6]-i[5])+' '+str(0)+')'
                str2.append(str1)
        elif savetype == 'circle':
            for i in data_new:
                str1 = 'circle('+str(i[1])+' '+str(i[2])+' '+str(i[3])+')'
                str2.append(str1)
        
        with open('a.reg','a+',encoding='utf-8') as f:
            for i in str2:
                f.write(i+'\n')
            f.close()

#---------------------------------------------------------------convolve6-----------------------------------------------------------------:
if choisemod == 6:
    def convolve(n,kn):
        hdu1 = fits.open(fitspath1+fitsname1)
        hdu2 = fits.open(fitspath2+fitsname2)
        
        img_change = hdu1[n].data
        hdr_change = hdu1[n].header
        img2 = hdu2[kn].data
        
        f = signal.fftconvolve(img_change,img2,mode='same')
        grey_change = fits.ImageHDU(f,hdr_change)
        
        if len(hdu1)>0:
            hdr10=hdu1[0].header
            img10=hdu1[0].data
            grey0=fits.PrimaryHDU(img10,hdr10)
            greylist=[grey0]
        if len(hdu1)>1:
            hdr11=hdu1[1].header
            img11=hdu1[1].data
            grey1=fits.ImageHDU(img11,hdr11)
            greylist=[grey0,grey1]
        if len(hdu1)>2:
            hdr12=hdu1[2].header
            img12=hdu1[2].data
            grey2=fits.ImageHDU(img12,hdr12)
            greylist=[grey0,grey1,grey2]
        if len(hdu1)>3:
            hdr13=hdu1[3].header
            img13=hdu1[3].data
            grey3=fits.ImageHDU(img13,hdr13)
            greylist=[grey0,grey1,grey2,grey3]
        
        if n == 0:
            greylist[0]=fits.PrimaryHDU(f,hdr_change)
        if n == 1:
            greylist[1]=grey_change
        if n == 2:
            greylist[2]=grey_change
        if n == 3:
            greylist[3]=grey_change
        greyHDU=fits.HDUList(greylist)
        greyHDU.writeto('newfits.fits')

#----------------------------------------------------------------usemodel-----------------------------------------------------------------:
if __name__ == '__main__':
    if choisemod == 0:
        kpctopix(z_peak,pixscale)
        print(pix_per1kpc)
    if choisemod == 0.1:
        savecsv(name,colname,data)
    if choisemod == 0.2:
        cutfits(hdu,mode,center,width)
    if choisemod == 0.3:
        savefits(data)
    if choisemod == 0.4:
        MJysr_to_Jypix(pixscal,SB)
    if choisemod == 1:
        pltfits1(Taxesset=Taxesset)
    if choisemod == 2:
        try:pltfits2(csvname1,csvname2=csvname2,Taxesset=Taxesset)
        except:pltfits2(csvname1,Taxesset=Taxesset)
    if choisemod == 3:
        pltfits3(csvname,xytype,numrows,numcols,picscale)
    if choisemod == 4:
        searchfilename()
        print(hudfp,xdf,hff)
    if choisemod == 5:
        pltfits5(csvpath,csvname,selectname,savetype)
    if choisemod == 6:
        convolve(n,kn)
