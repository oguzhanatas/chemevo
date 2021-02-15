# -*- coding: utf-8 -*-
#subat 2021

from pandas import read_csv
import numpy as np
import pylab as py
import os
import matplotlib.pyplot as plt
from astropy.io import fits as pf
import warnings
import pandas as pd



warnings.simplefilter('ignore', np.RankWarning)
"""
#ram problemi icin.......
#https://stackoverflow.com/questions/41105733/limit-ram-usage-to-python-program
import resource
import platform
import sys

def memory_limit(percentage: float):
  
    #linux
    
    if platform.system() != "Linux":
        print('Only works on linux!')
        return
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (get_memory() * 1024 * percentage, hard))

def get_memory():
    with open('/proc/meminfo', 'r') as mem:
        free_memory = 0
        for i in mem:
            sline = i.split()
            if str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
                free_memory += int(sline[1])
    return free_memory

def memory(percentage=0.8):
    def decorator(function):
        def wrapper(*args, **kwargs):
            memory_limit(percentage)
            try:
                function(*args, **kwargs)
            except MemoryError:
                mem = get_memory() / 1024 /1024
                print('Remain: %.2f GB' % mem)
                sys.stderr.write('\n\nERROR: Memory Exception\n')
                sys.exit(1)
        return wrapper
    return decorator

@memory(percentage=0.8)
def main():
    print('My memory is limited to 80%.')

"""





data=os.pardir+"/data/"
output=os.pardir+"/output/"
figures=os.pardir+"/figures/"

#apogee dr16 https://www.sdss.org/dr16/irspec/catalogs/ good sample getting
apogee_hdu = pf.open(data+"allStar-r12-l33.fits")
apogee_data =apogee_hdu[1].data
apogee_hdu.close()


starbad = 2**23 #bit flag for bad stars 
gd = np.bitwise_and(apogee_data["ASPCAPFLAG"], starbad) == 0 
teff_logg_check = np.logical_and(apogee_data["TEFF"] > 0, apogee_data["LOGG"] > -10) # this checks for -9999 values
teff_logg_feh_check = np.logical_and(apogee_data["FE_H"]> -6, teff_logg_check)

indices = np.where(np.logical_and(gd, teff_logg_feh_check)) 
apogee_good = apogee_data[indices] # this only the good data now


#galah DR3
gDR3_hdu = pf.open(data+'GALAH_DR3_main_allspec_v1.fits')
gDR3=gDR3_hdu[1].data
#cols_gDR3 = gDR3_hdu[1].columns
gDR3_hdu.close()

#print(cols_gDR3.names)



gd_galah = np.bitwise_and(gDR3["flag_sp"], 0) == 0
teff_logg_check_galah = np.logical_and(gDR3["teff"] > 0, gDR3["logg"] > -10) # this checks for -9999 values
teff_logg_feh_check_galah = np.logical_and(gDR3["fe_h"]> -6, teff_logg_check_galah)
indices_galah = np.where(np.logical_and(gd_galah, teff_logg_feh_check_galah)) 
galah_good = gDR3[indices] # this only the good data now

#apogee datas https://data.sdss.org/datamodel/files/APOGEE_ASPCAP/APRED_VERS/ASPCAP_VERS/allStar.html
J=np.array(apogee_good['J'])
J_ERR=np.array(apogee_good['J_ERR'])
H=np.array(apogee_good['H'])
H_ERR=np.array(apogee_good['H_ERR'])
K=np.array(apogee_good['K'])
K_ERR=np.array(apogee_good['K_ERR'])

RA=np.array(apogee_good['RA'])
DEC=np.array(apogee_good['DEC'])
GLON=np.array(apogee_good['GLON'])
GLAT=np.array(apogee_good['GLAT'])

SNR=np.array(apogee_good['SNR'])

VHELIO_AVG=np.array(apogee_good['VHELIO_AVG'])
VERR=np.array(apogee_good['VERR'])

TEFF=np.array(apogee_good['TEFF'])
TEFF_ERR=np.array(apogee_good['TEFF_ERR'])
LOGG=np.array(apogee_good['LOGG'])
LOGG_ERR=np.array(apogee_good['LOGG_ERR'])


ALPHA_M_a=np.array(apogee_good['ALPHA_M'])
ALPHA_M_ERR_a=np.array(apogee_good['ALPHA_M_ERR'])
M_H_a=np.array(apogee_good['M_H'])
M_H_ERR_a=np.array(apogee_good['M_H_ERR'])
C_FE=np.array(apogee_good['C_FE'])
O_FE=np.array(apogee_good['O_FE'])
NA_FE=np.array(apogee_good['NA_FE'])
MG_FE=np.array(apogee_good['MG_FE'])
AL_FE=np.array(apogee_good['AL_FE'])
SI_FE=np.array(apogee_good['SI_FE'])
K_FE=np.array(apogee_good['K_FE'])
CA_FE=np.array(apogee_good['CA_FE'])
TI_FE=np.array(apogee_good['TI_FE'])
TIII_FE=np.array(apogee_good['TIII_FE'])
V_FE=np.array(apogee_good['V_FE'])
CR_FE=np.array(apogee_good['CR_FE'])
MN_FE=np.array(apogee_good['MN_FE'])
FE_H=np.array(apogee_good['FE_H'])
CO_FE=np.array(apogee_good['CO_FE'])
NI_FE=np.array(apogee_good['NI_FE'])
CU_FE=np.array(apogee_good['CU_FE'])
RB_FE=np.array(apogee_good['RB_FE'])
CE_FE=np.array(apogee_good['CE_FE'])
ND_FE=np.array(apogee_good['ND_FE'])
C_FE_ERR=np.array(apogee_good['C_FE_ERR'])
O_FE_ERR=np.array(apogee_good['O_FE_ERR'])
NA_FE_ERR=np.array(apogee_good['NA_FE_ERR'])
MG_FE_ERR=np.array(apogee_good['MG_FE_ERR'])
AL_FE_ERR=np.array(apogee_good['AL_FE_ERR'])
SI_FE_ERR=np.array(apogee_good['SI_FE_ERR'])
K_FE_ERR=np.array(apogee_good['K_FE_ERR'])
CA_FE_ERR=np.array(apogee_good['CA_FE_ERR'])
TI_FE_ERR=np.array(apogee_good['TI_FE_ERR'])
TIII_FE_ERR=np.array(apogee_good['TIII_FE_ERR'])
V_FE_ERR=np.array(apogee_good['V_FE_ERR'])
CR_FE_ERR=np.array(apogee_good['CR_FE_ERR'])
MN_FE_ERR=np.array(apogee_good['MN_FE_ERR'])
FE_H_ERR=np.array(apogee_good['FE_H_ERR'])
CO_FE_ERR=np.array(apogee_good['CO_FE_ERR'])
NI_FE_ERR=np.array(apogee_good['NI_FE_ERR'])
CU_FE_ERR=np.array(apogee_good['CU_FE_ERR'])
RB_FE_ERR=np.array(apogee_good['RB_FE_ERR'])
CE_FE_ERR=np.array(apogee_good['CE_FE_ERR'])
ND_FE_ERR=np.array(apogee_good['ND_FE_ERR'])



GAIA_SOURCE_ID=np.array(apogee_good['GAIA_SOURCE_ID'])
GAIA_PARALLAX=np.array(apogee_good['GAIA_PARALLAX'])
GAIA_PARALLAX_ERROR=np.array(apogee_good['GAIA_PARALLAX_ERROR'])
GAIA_PMRA=np.array(apogee_good['GAIA_PMRA'])
GAIA_PMRA_ERROR=np.array(apogee_good['GAIA_PMRA_ERROR'])
GAIA_PMDEC=np.array(apogee_good['GAIA_PMDEC'])
GAIA_PMDEC_ERROR=np.array(apogee_good['GAIA_PMDEC_ERROR'])
GAIA_PHOT_G_MEAN_MAG=np.array(apogee_good['GAIA_PHOT_G_MEAN_MAG'])
GAIA_PHOT_BP_MEAN_MAG=np.array(apogee_good['GAIA_PHOT_BP_MEAN_MAG'])
GAIA_PHOT_RP_MEAN_MAG=np.array(apogee_good['GAIA_PHOT_RP_MEAN_MAG'])

GAIA_PHOT_BP_RP_MEAN_MAG1=GAIA_PHOT_BP_MEAN_MAG-GAIA_PHOT_RP_MEAN_MAG
GAIA_PHOT_BP_RP_MEAN_MAG=np.array(GAIA_PHOT_BP_RP_MEAN_MAG1)

GAIA_RADIAL_VELOCITY=np.array(apogee_good['GAIA_RADIAL_VELOCITY'])
GAIA_RADIAL_VELOCITY_ERROR=np.array(apogee_good['GAIA_RADIAL_VELOCITY_ERROR'])









#galah datas https://docs.datacentral.org.au/galah/dr3/table-schema/
j_m=np.array(galah_good['j_m'])
j_msigcom=np.array(galah_good['j_msigcom'])
h_m=np.array(galah_good['h_m'])
h_msigcom=np.array(galah_good['h_msigcom'])
ks_m=np.array(galah_good['ks_m'])
ks_msigcom=np.array(galah_good['ks_msigcom'])


ra=np.array(galah_good['ra'])
dec=np.array(galah_good['dec'])
l=np.array(galah_good['l'])
b=np.array(galah_good['b'])


snr_c1=np.array(galah_good['snr_c1_iraf'])
snr_c2=np.array(galah_good['snr_c2_iraf'])
snr_c3=np.array(galah_good['snr_c3_iraf'])
snr_c4=np.array(galah_good['snr_c4_iraf'])

rv_galah=np.array(galah_good['rv_galah'])
e_rv_galah=np.array(galah_good['e_rv_galah'])

teff=np.array(galah_good['teff'])
e_teff=np.array(galah_good['e_teff'])
logg=np.array(galah_good['logg'])
e_logg=np.array(galah_good['e_logg'])

alpha_fe=np.array(galah_good['alpha_fe'])
e_alpha_fe=np.array(galah_good['e_alpha_fe'])



C_fe=np.array(galah_good['C_fe'])
e_C_fe=np.array(galah_good['e_C_fe'])
O_fe=np.array(galah_good['O_fe'])
e_O_fe=np.array(galah_good['e_O_fe'])
Na_fe=np.array(galah_good['Na_fe'])
e_Na_fe=np.array(galah_good['e_Na_fe'])
Mg_fe=np.array(galah_good['Mg_fe'])
e_Mg_fe=np.array(galah_good['e_Mg_fe'])
Al_fe=np.array(galah_good['Al_fe'])
e_Al_fe=np.array(galah_good['e_Al_fe'])
Si_fe=np.array(galah_good['Si_fe'])
e_Si_fe=np.array(galah_good['e_Si_fe'])
K_fe=np.array(galah_good['K_fe'])
e_K_fe=np.array(galah_good['e_K_fe'])
Ca_fe=np.array(galah_good['Ca_fe'])
e_Ca_fe=np.array(galah_good['e_Ca_fe'])
Ti_fe=np.array(galah_good['Ti_fe'])
e_Ti_fe=np.array(galah_good['e_Ti_fe'])
Ti2_fe=np.array(galah_good['Ti2_fe'])
e_Ti2_fe=np.array(galah_good['e_Ti2_fe'])
V_fe=np.array(galah_good['V_fe'])
e_V_fe=np.array(galah_good['e_V_fe'])
Cr_fe=np.array(galah_good['Cr_fe'])
e_Cr_fe=np.array(galah_good['e_Cr_fe'])
Mn_fe=np.array(galah_good['Mn_fe'])
e_Mn_fe=np.array(galah_good['e_Mn_fe'])
fe_h=np.array(galah_good['fe_h'])
e_fe_h=np.array(galah_good['e_fe_h'])
Co_fe=np.array(galah_good['Co_fe'])
e_Co_fe=np.array(galah_good['e_Co_fe'])
Ni_fe=np.array(galah_good['Ni_fe'])
e_Ni_fe=np.array(galah_good['e_Ni_fe'])
Cu_fe=np.array(galah_good['Cu_fe'])
e_Cu_fe=np.array(galah_good['e_Cu_fe'])
Rb_fe=np.array(galah_good['Rb_fe'])
e_Rb_fe=np.array(galah_good['e_Rb_fe'])
Ce_fe=np.array(galah_good['Ce_fe'])
e_Ce_fe=np.array(galah_good['e_Ce_fe'])
Nd_fe=np.array(galah_good['Nd_fe'])
e_Nd_fe=np.array(galah_good['e_Nd_fe'])

source_id=np.array(galah_good['source_id'])
parallax=np.array(galah_good['parallax'])
parallax_error=np.array(galah_good['parallax_error'])
pmra=np.array(galah_good['pmra'])
pmra_error=np.array(galah_good['pmra_error'])
pmdec=np.array(galah_good['pmdec'])
pmdec_error=np.array(galah_good['pmdec_error'])
phot_g_mean_mag=np.array(galah_good['phot_g_mean_mag'])
bp_rp=np.array(galah_good['bp_rp'])
rv_gaia=np.array(galah_good['rv_gaia'])
e_rv_gaia=np.array(galah_good['e_rv_gaia'])



#merge data............
j_all=np.concatenate((J,j_m),axis=0)
e_j_all=np.concatenate((J_ERR,j_msigcom),axis=0)
h_all=np.concatenate((H,h_m),axis=0)
e_h_all=np.concatenate((H_ERR,h_msigcom),axis=0)
k_all=np.concatenate((K,ks_m),axis=0)
e_k_all=np.concatenate((K_ERR,ks_msigcom),axis=0)

ra_all=np.concatenate((RA,ra),axis=0)
dec_all=np.concatenate((DEC,dec),axis=0)
l_all=np.concatenate((GLON,l),axis=0)
b_all=np.concatenate((GLAT,b),axis=0)

snr_a=np.concatenate((SNR,[-9999.0]*len(snr_c1)),axis=0)
snr_c1=np.concatenate(([-9999.0]*len(SNR),snr_c1),axis=0)
snr_c2=np.concatenate(([-9999.0]*len(SNR),snr_c2),axis=0)
snr_c3=np.concatenate(([-9999.0]*len(SNR),snr_c3),axis=0)
snr_c4=np.concatenate(([-9999.0]*len(SNR),snr_c4),axis=0)

rv_all=np.concatenate((VHELIO_AVG,rv_galah),axis=0)
e_rv_all=np.concatenate((VERR,e_rv_galah),axis=0)

teff_all=np.concatenate((TEFF,teff),axis=0)
e_teff_all=np.concatenate((TEFF_ERR,e_teff),axis=0)
logg_all=np.concatenate((LOGG,logg),axis=0)
e_logg_all=np.concatenate((LOGG_ERR,e_logg),axis=0)

ALPHA_M_a=np.concatenate((ALPHA_M_a,[-9999.0]*len(alpha_fe)),axis=0)
ALPHA_M_ERR_a=np.concatenate((ALPHA_M_ERR_a,[-9999.0]*len(alpha_fe)),axis=0)
M_H_a=np.concatenate((M_H_a,[-9999.0]*len(alpha_fe)),axis=0)
M_H_ERR_a=np.concatenate((M_H_ERR_a,[-9999.0]*len(alpha_fe)),axis=0)

alphafe_g=np.concatenate(([-9999.0]*len(ALPHA_M_a),alpha_fe),axis=0)
e_alphafe_g=np.concatenate(([-9999.0]*len(ALPHA_M_a),e_alpha_fe),axis=0)

cfe_all=np.concatenate((C_FE,C_fe),axis=0)
e_cfe_all=np.concatenate((C_FE_ERR,e_C_fe),axis=0)
ofe_all=np.concatenate((O_FE,O_fe),axis=0)
e_ofe_all=np.concatenate((O_FE_ERR,e_O_fe),axis=0)
nafe_all=np.concatenate((NA_FE,Na_fe),axis=0)
e_nafe_all=np.concatenate((NA_FE_ERR,e_Na_fe),axis=0)
mgfe_all=np.concatenate((MG_FE,Mg_fe),axis=0)
e_mgfe_all=np.concatenate((MG_FE_ERR,e_Mg_fe),axis=0)
alfe_all=np.concatenate((AL_FE,Al_fe),axis=0)
e_alfe_all=np.concatenate((AL_FE_ERR,e_Al_fe),axis=0)
sife_all=np.concatenate((SI_FE,Si_fe),axis=0)
e_sife_all=np.concatenate((SI_FE_ERR,e_Si_fe),axis=0)
kfe_all=np.concatenate((K_FE,K_fe),axis=0)
e_kfe_all=np.concatenate((K_FE_ERR,e_K_fe),axis=0)
cafe_all=np.concatenate((CA_FE,Ca_fe),axis=0)
e_cafe_all=np.concatenate((CA_FE_ERR,e_Ca_fe),axis=0)
tife_all=np.concatenate((TI_FE,Ti_fe),axis=0)
e_tife_all=np.concatenate((TI_FE_ERR,e_Ti_fe),axis=0)
ti2fe_all=np.concatenate((TIII_FE,Ti2_fe),axis=0)
e_ti2fe_all=np.concatenate((TIII_FE_ERR,e_Ti2_fe),axis=0)
vfe_all=np.concatenate((V_FE,V_fe),axis=0)
e_vfe_all=np.concatenate((V_FE_ERR,e_V_fe),axis=0)
crfe_all=np.concatenate((CR_FE,Cr_fe),axis=0)
e_crfe_all=np.concatenate((CR_FE_ERR,e_Cr_fe),axis=0)
mnfe_all=np.concatenate((MN_FE,Mn_fe),axis=0)
e_mnfe_all=np.concatenate((MN_FE_ERR,e_Mn_fe),axis=0)
feh_all=np.concatenate((FE_H,fe_h),axis=0)
e_feh_all=np.concatenate((FE_H_ERR,e_fe_h),axis=0)
cofe_all=np.concatenate((CO_FE,Co_fe),axis=0)
e_cofe_all=np.concatenate((CO_FE_ERR,e_Co_fe),axis=0)
nife_all=np.concatenate((NI_FE,Ni_fe),axis=0)
e_nife_all=np.concatenate((NI_FE_ERR,e_Ni_fe),axis=0)
cufe_all=np.concatenate((CU_FE,Cu_fe),axis=0)
e_cufe_all=np.concatenate((CU_FE_ERR,e_Cu_fe),axis=0)
rbfe_all=np.concatenate((RB_FE,Rb_fe),axis=0)
e_rbfe_all=np.concatenate((RB_FE_ERR,e_Rb_fe),axis=0)
cefe_all=np.concatenate((CE_FE,Ce_fe),axis=0)
e_cefe_all=np.concatenate((CE_FE_ERR,e_Ce_fe),axis=0)
ndfe_all=np.concatenate((ND_FE,Nd_fe),axis=0)
e_ndfe_all=np.concatenate((ND_FE_ERR,e_Nd_fe),axis=0)

gaia_soruceid_all=np.concatenate((GAIA_SOURCE_ID,source_id),axis=0)
gaia_parallax_all=np.concatenate((GAIA_PARALLAX,parallax),axis=0)
e_gaia_parallax_all=np.concatenate((GAIA_PARALLAX_ERROR,parallax_error),axis=0)
gaia_pmra_all=np.concatenate((GAIA_PMRA,pmra),axis=0)
e_gaia_pmra_all=np.concatenate((GAIA_PMRA_ERROR,pmra_error),axis=0)
gaia_pmdec_all=np.concatenate((GAIA_PMDEC,pmdec),axis=0)
e_gaia_pmdec_all=np.concatenate((GAIA_PMDEC_ERROR,pmdec_error),axis=0)
gaia_gmag_all=np.concatenate((GAIA_PHOT_G_MEAN_MAG,phot_g_mean_mag),axis=0)
gaia_bprp_all=np.concatenate((GAIA_PHOT_BP_RP_MEAN_MAG,bp_rp),axis=0)
gaia_rv_all=np.concatenate((GAIA_RADIAL_VELOCITY,rv_gaia),axis=0)
e_gaia_rv_all=np.concatenate((GAIA_RADIAL_VELOCITY_ERROR,e_rv_gaia),axis=0)

f=open(output+"apogee_galah_good.csv","w")

f.writelines('j_all,e_j_all,h_all,e_h_all,k_all,e_k_all,ra_all,dec_all,l_all,b_all,snr_a,snr_c1,snr_c2,snr_c3,snr_c4,rv_all,e_rv_all,teff_all,e_teff_all,logg_all,e_logg_all,ALPHA_M_a,ALPHA_M_ERR_a,M_H_a,M_H_ERR_a,alphafe_g,e_alphafe_g,cfe_all,e_cfe_all,ofe_all,e_ofe_all,nafe_all,e_nafe_all,mgfe_all,e_mgfe_all,alfe_all,e_alfe_all,sife_all,e_sife_all,kfe_all,e_kfe_all,cafe_all,e_cafe_all,tife_all,e_tife_all,ti2fe_all,e_ti2fe_all,vfe_all,e_vfe_all,crfe_all,e_crfe_all,mnfe_all,e_mnfe_all,feh_all,e_feh_all,cofe_all,e_cofe_all,nife_all,e_nife_all,cufe_all,e_cufe_all,rbfe_all,e_rbfe_all,cefe_all,e_cefe_all,ndfe_all,e_ndfe_all,gaia_soruceid_all,gaia_parallax_all,e_gaia_parallax_all,gaia_pmra_all,e_gaia_pmra_all,gaia_pmdec_all,e_gaia_pmdec_all,gaia_gmag_all,gaia_bprp_all,gaia_rv_all,e_gaia_rv_all'+'\n')
             
for i in range(len(j_all)):
	f.writelines(str(j_all[i])+','+str(e_j_all[i])+','+str(h_all[i])+','+str(e_h_all[i])+','+str(k_all[i])+','+str(e_k_all[i])+','+str(ra_all[i])+','+str(dec_all[i])+','+str(l_all[i])+','+str(b_all[i])+','+str(snr_a[i])+','+str(snr_c1[i])+','+str(snr_c2[i])+','+str(snr_c3[i])+','+str(snr_c4[i])+','+str(rv_all[i])+','+str(e_rv_all[i])+','+str(teff_all[i])+','+str(e_teff_all[i])+','+str(logg_all[i])+','+str(e_logg_all[i])+','+str(ALPHA_M_a[i])+','+str(ALPHA_M_ERR_a[i])+','+str(M_H_a[i])+','+str(M_H_ERR_a[i])+','+str(alphafe_g[i])+','+str(e_alphafe_g[i])+','+str(cfe_all[i])+','+str(e_cfe_all[i])+','+str(ofe_all[i])+','+str(e_ofe_all[i])+','+str(nafe_all[i])+','+str(e_nafe_all[i])+','+str(mgfe_all[i])+','+str(e_mgfe_all[i])+','+str(alfe_all[i])+','+str(e_alfe_all[i])+','+str(sife_all[i])+','+str(e_sife_all[i])+','+str(kfe_all[i])+','+str(e_kfe_all[i])+','+str(cafe_all[i])+','+str(e_cafe_all[i])+','+str(tife_all[i])+','+str(e_tife_all[i])+','+str(ti2fe_all[i])+','+str(e_ti2fe_all[i])+','+str(vfe_all[i])+','+str(e_vfe_all[i])+','+str(crfe_all[i])+','+str(e_crfe_all[i])+','+str(mnfe_all[i])+','+str(e_mnfe_all[i])+','+str(feh_all[i])+','+str(e_feh_all[i])+','+str(cofe_all[i])+','+str(e_cofe_all[i])+','+str(nife_all[i])+','+str(e_nife_all[i])+','+str(cufe_all[i])+','+str(e_cufe_all[i])+','+str(rbfe_all[i])+','+str(e_rbfe_all[i])+','+str(cefe_all[i])+','+str(e_cefe_all[i])+','+str(ndfe_all[i])+','+str(e_ndfe_all[i])+','+str(gaia_soruceid_all[i])+','+str(gaia_parallax_all[i])+','+str(e_gaia_parallax_all[i])+','+str(gaia_pmra_all[i])+','+str(e_gaia_pmra_all[i])+','+str(gaia_pmdec_all[i])+','+str(e_gaia_pmdec_all[i])+','+str(gaia_gmag_all[i])+','+str(gaia_bprp_all[i])+','+str(gaia_rv_all[i])+','+str(e_gaia_rv_all[i])+'\n')
f.close()
#print(max(teff_all),min(teff_all))
#print(max(logg_all),min(logg_all))
#print(max(feh_all),min(feh_all))
#print(max(alphafe_g),min(alphafe_g))
#print(max(ALPHA_M_a),min(ALPHA_M_a))
"""
#fits dosyasi yaratma.............
#https://docs.astropy.org/en/stable/io/fits/
#burada creating new table file başlığını oku....

from astropy.io import fits

c1=fits.Column(name='j_all',array=j_all,format='K')
c2=fits.Column(name='e_j_all',array=e_j_all,format='K')
c3=fits.Column(name='h_all',array=h_all,format='K')
c4=fits.Column(name='e_h_all',array=e_h_all,format='K')
c5=fits.Column(name='k_all',array=k_all,format='K')
c6=fits.Column(name='e_k_all',array=e_k_all,format='K')
c7=fits.Column(name='ra_all',array=ra_all,format='K')
c8=fits.Column(name='dec_all',array=dec_all,format='K')
c9=fits.Column(name='l_all',array=l_all,format='K')
c10=fits.Column(name='b_all',array=b_all,format='K')
c11=fits.Column(name='snr_a',array=snr_a,format='K')
c12=fits.Column(name='snr_c1',array=snr_c1,format='K')
c13=fits.Column(name='snr_c2',array=snr_c2,format='K')
c14=fits.Column(name='snr_c3',array=snr_c3,format='K')
c15=fits.Column(name='snr_c4',array=snr_c4,format='K')
c16=fits.Column(name='rv_all',array=rv_all,format='K')
c17=fits.Column(name='e_rv_all',array=e_rv_all,format='K')
c18=fits.Column(name='teff_all',array=teff_all,format='K')
c19=fits.Column(name='e_teff_all',array=e_teff_all,format='K')
c20=fits.Column(name='logg_all',array=logg_all,format='K')
c21=fits.Column(name='e_logg_all',array=e_logg_all,format='K')
#c22=fits.Column(name='alphafe_all',array=alphafe_all,format='K')
#c23=fits.Column(name='e_alphafe_all',array=e_alphafe_all,format='K')


c22_1=fits.Column(name="ALPHA_M_a",array=ALPHA_M_a,format='K')
c22_2=fits.Column(name="ALPHA_M_ERR_a",array=ALPHA_M_ERR_a,format='K')
c22_3=fits.Column(name="M_H_a",array=M_H_a,format='K')
c22_4=fits.Column(name="M_H_ERR_a",array=M_H_ERR_a,format='K')

c22_5=fits.Column(name="alphafe_g",array=alphafe_g,format='K')
c22_6=fits.Column(name="e_alphafe_g",array=e_alphafe_g,format='K')


c24=fits.Column(name='cfe_all',array=cfe_all,format='K')
c25=fits.Column(name='e_cfe_all',array=e_cfe_all,format='K')
c26=fits.Column(name='ofe_all',array=ofe_all,format='K')
c27=fits.Column(name='e_ofe_all',array=e_ofe_all,format='K')
c28=fits.Column(name='nafe_all',array=nafe_all,format='K')
c29=fits.Column(name='e_nafe_all',array=e_nafe_all,format='K')
c30=fits.Column(name='mgfe_all',array=mgfe_all,format='K')
c31=fits.Column(name='e_mgfe_all',array=e_mgfe_all,format='K')
c32=fits.Column(name='alfe_all',array=alfe_all,format='K')
c33=fits.Column(name='e_alfe_all',array=e_alfe_all,format='K')
c34=fits.Column(name='sife_all',array=sife_all,format='K')
c35=fits.Column(name='e_sife_all',array=e_sife_all,format='K')
c36=fits.Column(name='kfe_all',array=kfe_all,format='K')
c37=fits.Column(name='e_kfe_all',array=e_kfe_all,format='K')
c38=fits.Column(name='cafe_all',array=cafe_all,format='K')
c39=fits.Column(name='e_cafe_all',array=e_cafe_all,format='K')
c40=fits.Column(name='tife_all',array=tife_all,format='K')
c41=fits.Column(name='e_tife_all',array=e_tife_all,format='K')
c42=fits.Column(name='ti2fe_all',array=ti2fe_all,format='K')
c43=fits.Column(name='e_ti2fe_all',array=e_ti2fe_all,format='K')
c44=fits.Column(name='vfe_all',array=vfe_all,format='K')
c45=fits.Column(name='e_vfe_all',array=e_vfe_all,format='K')
c46=fits.Column(name='crfe_all',array=crfe_all,format='K')
c47=fits.Column(name='e_crfe_all',array=e_crfe_all,format='K')
c48=fits.Column(name='mnfe_all',array=mnfe_all,format='K')
c49=fits.Column(name='e_mnfe_all',array=e_mnfe_all,format='K')
c50=fits.Column(name='feh_all',array=feh_all,format='K')
c51=fits.Column(name='e_feh_all',array=e_feh_all,format='K')
c52=fits.Column(name='cofe_all',array=cofe_all,format='K')
c53=fits.Column(name='e_cofe_all',array=e_cofe_all,format='K')
c54=fits.Column(name='nife_all',array=nife_all,format='K')
c55=fits.Column(name='e_nife_all',array=e_nife_all,format='K')
c56=fits.Column(name='cufe_all',array=cufe_all,format='K')
c57=fits.Column(name='e_cufe_all',array=e_cufe_all,format='K')
c58=fits.Column(name='rbfe_all',array=rbfe_all,format='K')
c59=fits.Column(name='e_rbfe_all',array=e_rbfe_all,format='K')
c60=fits.Column(name='cefe_all',array=cefe_all,format='K')
c61=fits.Column(name='e_cefe_all',array=e_cefe_all,format='K')
c62=fits.Column(name='ndfe_all',array=ndfe_all,format='K')
c63=fits.Column(name='e_ndfe_all',array=e_ndfe_all,format='K')
c64=fits.Column(name='g_si',array=gaia_soruceid_all,format='K')
c65=fits.Column(name='g_pl',array=gaia_parallax_all,format='K')
c66=fits.Column(name='e_g_pl',array=e_gaia_parallax_all,format='K')
c67=fits.Column(name='g_pmra',array=gaia_pmra_all,format='K')
c68=fits.Column(name='e_g_pmra',array=e_gaia_pmra_all,format='K')
c69=fits.Column(name='g_pmdec',array=gaia_pmdec_all,format='K')
c70=fits.Column(name='e_g_pmdec',array=e_gaia_pmdec_all,format='K')
c71=fits.Column(name='g_gmag',array=gaia_gmag_all,format='K')
c72=fits.Column(name='g_bprp',array=gaia_bprp_all,format='K')
c73=fits.Column(name='g_rv',array=gaia_rv_all,format='K')
c74=fits.Column(name='e_g_rv',array=e_gaia_rv_all,format='K')

t=fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22_1,c22_2,c22_3,c22_4,c22_5,c22_6,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,c38,c39,c40,c41,c42,c43,c44,c45,c46,c47,c48,c49,c50,c51,c52,c53,c54,c55,c56,c57,c58,c59,c60,c61,c62,c63,c64,c65,c66,c67,c68,c69,c70,c71,c72,c73,c74])
t.writeto(output+'apogee_galah_good.fits')

"""


