from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
import pandas as pd

startup_time=[]
hour=[]
minute=[]
target_list_path = 'speculoos_target_list_v6.txt'

def Path_txt_files(telescope):
    Path=os.path.join('./DATABASE',telescope,'Plans_by_date')
    return Path

def startup_no_flats(t_now,name,sun_set,date_start,Path):
    date_start=np.datetime64(date_start)
    sun_set=np.datetime64(sun_set)
    hour00=sun_set.astype(object).hour
    hour00='{:02d}'.format(int(hour00))
    minute00=sun_set.astype(object).minute
    minute00='{:02d}'.format(int(minute00))
    hour0=date_start.astype(object).hour
    hour0='{:02d}'.format(int(hour0))
    minute0=date_start.astype(object).minute
    minute0='{:02d}'.format(int(minute0))
    startup_time1=sun_set+np.timedelta64(10,'m')
    hour1=startup_time1.astype(object).hour
    hour1='{:02d}'.format(int(hour1))
    minute1=startup_time1.astype(object).minute
    minute1='{:02d}'.format(int(minute1))
    startup_time2=date_start-np.timedelta64(40,'m')
    hour2=startup_time2.astype(object).hour
    hour2='{:02d}'.format(int(hour2))
    minute2=startup_time2.astype(object).minute
    minute2='{:02d}'.format(int(minute2))
    startup_time3=date_start-np.timedelta64(50,'m')
    hour3=startup_time3.astype(object).hour
    hour3='{:02d}'.format(int(hour3))
    minute3=startup_time3.astype(object).minute
    minute3='{:02d}'.format(int(minute3))
    with open(os.path.join(Path,str(t_now),'startup.txt'),'w') as out:
        str00=';\n'
        str0=' == Startup =='
        str1='#waituntil 1, '
        str2='#domeopen\n'
        str3='#chill '
        str33='-60\n'
        str4=';#duskflats Cal_flatexo.txt\n'
        str5='#chain '
        str6='#quitat '
        out.write(';' + str0 + '\n')
        for i in range(1,3):
            out.write(str00)
        out.write(str1 + pd.to_datetime(str(sun_set)).strftime('%Y/%m/%d %H:%M:%S') + '\n')#str(hour00) + ':' + str(minute00) + '\n')
        out.write(str2)
        for i in range(1,2):
            out.write(str00)
        #out.write(str1 + str(hour2) + ':' + str(minute2) + '\n')
        out.write(str3 + str33)
        for i in range(1,2):
            out.write(str00)
        out.write(str1 + pd.to_datetime(str(startup_time1)).strftime('%Y/%m/%d %H:%M:%S') + '\n')#str(hour1) + ':' + str(minute1) + '\n')
        out.write(str4)
        for i in range(1,2):
            out.write(str00)
        out.write(str6 + pd.to_datetime(str(date_start)).strftime('%Y/%m/%d %H:%M:%S') + '\n')#str(hour0) + ':' + str(minute0) + '\n')
        for i in range(1,2):
            out.write(str00)
        out.write(str5 + 'Obj_' + name + '.txt' + '\n')
        out.write(str00)

def startup(t_now,name,sun_set,date_start,Path):
    date_start=np.datetime64(date_start)
    sun_set=np.datetime64(sun_set)
    hour00=sun_set.astype(object).hour
    hour00='{:02d}'.format(int(hour00))
    minute00=sun_set.astype(object).minute
    minute00='{:02d}'.format(int(minute00))
    hour0=date_start.astype(object).hour
    hour0='{:02d}'.format(int(hour0))
    minute0=date_start.astype(object).minute
    minute0='{:02d}'.format(int(minute0))
    startup_time1=sun_set+np.timedelta64(10,'m')
    hour1=startup_time1.astype(object).hour
    hour1='{:02d}'.format(int(hour1))
    minute1=startup_time1.astype(object).minute
    minute1='{:02d}'.format(int(minute1))
    startup_time2=date_start-np.timedelta64(40,'m')
    hour2=startup_time2.astype(object).hour
    hour2='{:02d}'.format(int(hour2))
    minute2=startup_time2.astype(object).minute
    minute2='{:02d}'.format(int(minute2))
    startup_time3=date_start-np.timedelta64(50,'m')
    hour3=startup_time3.astype(object).hour
    hour3='{:02d}'.format(int(hour3))
    minute3=startup_time3.astype(object).minute
    minute3='{:02d}'.format(int(minute3))
    with open(os.path.join(Path,str(t_now),'startup.txt'),'w') as out:
        str00=';\n'
        str0=' == Startup =='
        str1='#waituntil 1, '
        str2='#domeopen\n'
        str3='#chill '
        str33='-60\n'
        str4='#duskflats Cal_flatexo.txt\n'
        str5='#chain '
        str6='#quitat '
        out.write(';' + str0 + '\n')
        for i in range(1,3):
            out.write(str00)
        out.write(str1 + pd.to_datetime(str(sun_set)).strftime('%Y/%m/%d %H:%M:%S') + '\n')#+ str(hour00) + ':' + str(minute00) + '\n')
        out.write(str2)
        for i in range(1,2):
            out.write(str00)
        #out.write(str1 + str(hour2) + ':' + str(minute2) + '\n')
        out.write(str3 + str33)
        for i in range(1,2):
            out.write(str00)
        out.write(str1 + pd.to_datetime(str(startup_time1)).strftime('%Y/%m/%d %H:%M:%S') + '\n')#str(hour1) + ':' + str(minute1) + '\n')
        out.write(str4)
        for i in range(1,2):
            out.write(str00)
        out.write(str6 + pd.to_datetime(str(date_start)).strftime('%Y/%m/%d %H:%M:%S') + '\n')#str(hour0) + ':' + str(minute0) + '\n')
        for i in range(1,2):
            out.write(str00)
        out.write(str5 + 'Obj_' + name + '.txt' + '\n')
        out.write(str00)

def startup_artemis(t_now,name,sun_set,date_start,Path):
    date_start=np.datetime64(date_start)
    sun_set=np.datetime64(sun_set)
    hour00=sun_set.astype(object).hour
    hour00='{:02d}'.format(int(hour00))
    minute00=sun_set.astype(object).minute
    minute00='{:02d}'.format(int(minute00))
    hour0=date_start.astype(object).hour
    hour0='{:02d}'.format(int(hour0))
    minute0=date_start.astype(object).minute
    minute0='{:02d}'.format(int(minute0))
    startup_time1=sun_set+np.timedelta64(10,'m')
    hour1=startup_time1.astype(object).hour
    hour1='{:02d}'.format(int(hour1))
    minute1=startup_time1.astype(object).minute
    minute1='{:02d}'.format(int(minute1))
    startup_time2=date_start-np.timedelta64(40,'m')
    hour2=startup_time2.astype(object).hour
    hour2='{:02d}'.format(int(hour2))
    minute2=startup_time2.astype(object).minute
    minute2='{:02d}'.format(int(minute2))
    startup_time3=date_start-np.timedelta64(50,'m')
    hour3=startup_time3.astype(object).hour
    hour3='{:02d}'.format(int(hour3))
    minute3=startup_time3.astype(object).minute
    minute3='{:02d}'.format(int(minute3))
    with open(os.path.join(Path,str(t_now),'startup.txt'),'w') as out:
        str00=';\n'
        str0=' == Startup =='
        str1='#waituntil 1, '
        str2='#domeopen\n'
        str3='#chill '
        str33='-55\n'
        str4='#duskflats Cal_flatexo_evening.txt\n'
        str5='#chain '
        str6='#quitat '
        out.write(';' + str0 + '\n')
        for i in range(1,3):
            out.write(str00)
        out.write(str1 + pd.to_datetime(str(sun_set)).strftime('%Y/%m/%d %H:%M:%S') + '\n')#+ str(hour00) + ':' + str(minute00) + '\n')
        out.write(str2)
        for i in range(1,2):
            out.write(str00)
        #out.write(str1 + str(hour2) + ':' + str(minute2) + '\n')
        out.write(str3 + str33)
        for i in range(1,2):
            out.write(str00)
        out.write(str1 + pd.to_datetime(str(startup_time1)).strftime('%Y/%m/%d %H:%M:%S') + '\n')#str(hour1) + ':' + str(minute1) + '\n')
        out.write(str4)
        for i in range(1,2):
            out.write(str00)
        out.write(str6 + pd.to_datetime(str(date_start)).strftime('%Y/%m/%d %H:%M:%S') + '\n')#str(hour0) + ':' + str(minute0) + '\n')
        for i in range(1,2):
            out.write(str00)
        out.write(str5 + 'Obj_' + name + '.txt' + '\n')
        out.write(str00)

def target(t_now,name,date_start,date_end,waitlimit,afinterval,autofocus,count,filt,exptime,ra1,ra2,ra3,dec1,dec2,dec3,name_2,Path,telescope=None):
    df = pd.read_csv(target_list_path,delimiter = ' ',index_col = False)
    #print(name)
    idx_target = None
    if 'Trappist' in name:
        idx_target = np.where((df['Sp_ID'] == 'Sp2306-0502'))[0]
        gaia_id_target = int(df['Gaia_ID'][idx_target].values)
    elif 'Sp' in name:
        df2 = pd.read_csv(target_list_path,delimiter=' ',index_col=None)
        #idx_target = np.where((df['spc'] == name.replace('_2','')))[0]
        idx_target = np.where((df2['Sp_ID'] == name.replace('_2','')))[0]
        gaia_id_target = df2['Gaia_ID'][int(idx_target)]
    if idx_target is None:
        df2 = pd.read_csv('target_list_special.txt',delimiter=' ',index_col=None)
        idx_target = np.where((df2['Sp_ID'] == name))[0]
        try:
            gaia_id_target = df2['Gaia_ID'][int(idx_target)]
        except TypeError:
            gaia_id_target = 'None'
    date_start=np.datetime64(date_start)
    date_end=np.datetime64(date_end)
    binning=1
    name_file=name + '.txt'
    hour0=date_start.astype(object).hour
    hour0='{:02d}'.format(int(hour0))
    minute0=date_start.astype(object).minute
    minute0='{:02d}'.format(int(minute0))
    hour1=date_end.astype(object).hour
    hour1='{:02d}'.format(int(hour1))
    minute1=date_end.astype(object).minute
    minute1='{:02d}'.format(int(minute1))
    with open(os.path.join(Path,str(t_now),'Obj_' + name_file),'w') as out:
        str00=';'
        str1='#waituntil 1, '
        str2='#waitinlimits '
        str22='#domeopen\n'
        str3='#nopreview \n'
        str33=';#afinterval '
        str4='#autofocus \n'
        str5='#count '
        str6='#binning '
        str7='#filter '
        str8='#interval '
        str9='#quitat '
        str10='#chain '
        s=''
        s2=''
        seq_ra=(str(int(ra1)),'h',str(int(ra2)),'m',str(ra3),'s')
        seq_dec=(str(int(dec1)),'d',str(abs(int(dec2))),'m',str(abs(dec3)),'s')
        s.join( seq_ra )
        s2.join( seq_dec )
        #print(s.join( seq_ra ),s2.join( seq_dec ))
        c = SkyCoord(s.join( seq_ra ),s2.join( seq_dec ),frame='icrs')
        #print(c.dec.degree,c.dec.radian)
        for i in range(1,2):
            out.write(str00 + '\n')
        out.write(str00 + ' ' + str(name).replace('_2','') + '\n')
        out.write(str00 + '\n')
        out.write(str00 + ' ' + str(gaia_id_target) + '\n')
        for i in range(1,2):
            out.write(str00 + '\n')
        out.write(str1 + str(hour0) + ':' + str(minute0) + '\n')
        out.write(str22)
        if telescope == 'Artemis':
            out.write('#chill -55\n')
        else:
            out.write('#chill -60\n')
        out.write(str2 + str(waitlimit) + '\n')
        out.write(str3)
        out.write(str33 + str(afinterval) +'\n')
        if autofocus is None:
            out.write(str00 + str4)
        out.write(str5 + str(count) + '\n')
        out.write(str6 + str(binning) + '\n')
        out.write(str7 + str(filt) + '\n')
        out.write(str8 + str(exptime) +'\n')
        out.write('#TAG Donuts=on'+'\n')
        if int(dec1)==-0:
            #print('ici')
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n') #8*np.cos(c.dec.radian)
        if int(dec1)<0 and int(dec1)!=-0:
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(abs(dec1)))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        if int(dec1)>0:
            #print('la')
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '+' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        out.write(str00 + '\n')
        out.write(str9 + str(hour1) + ':' + str(minute1) + '\n' )
        if name_2 is None:
            out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
        else:
            out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
        out.write(str00 + '\n')

def target_no_DONUTS(t_now,name,date_start,date_end,waitlimit,afinterval,autofocus,count,filt,exptime,ra1,ra2,ra3,dec1,dec2,dec3,name_2,Path):
    df = pd.read_csv(target_list_path,delimiter = ' ',index_col = False)
    idx_target = None
    if 'Trappist' in name:
        idx_target = np.where((df['Sp_ID'] == 'Sp2306-0502'))[0]
        gaia_id_target = int(df['Gaia_ID'][idx_target].values)
    elif 'Sp' in name:
        df2 = pd.read_csv(target_list_path,delimiter=' ',index_col=None)
        #idx_target = np.where((df['spc'] == name.replace('_2','')))[0]
        idx_target = np.where((df2['Sp_ID'] == name.replace('_2','')))[0]
        gaia_id_target = df2['Gaia_ID'][int(idx_target)] #int(df['gaia'][idx_target].values)
    if idx_target is None:
        df2 = pd.read_csv('target_list_special.txt',delimiter=' ',index_col=None)
        idx_target = np.where((df2['Sp_ID'] == name))[0]
        gaia_id_target = df2['Gaia_ID'][int(idx_target)]

    date_start=np.datetime64(date_start)
    date_end=np.datetime64(date_end)
    binning=1
    name_file=name + '.txt'
    hour0=date_start.astype(object).hour
    hour0='{:02d}'.format(int(hour0))
    minute0=date_start.astype(object).minute
    minute0='{:02d}'.format(int(minute0))
    hour1=date_end.astype(object).hour
    hour1='{:02d}'.format(int(hour1))
    minute1=date_end.astype(object).minute
    minute1='{:02d}'.format(int(minute1))
    with open(os.path.join(Path,str(t_now),'Obj_' + name_file),'w') as out:
        str00=';'
        str1='#waituntil 1, '
        str2='#waitinlimits '
        str22='#domeopen\n'
        str3='#nopreview \n'
        str33=';#afinterval '
        str4='#autofocus \n'
        str5='#count '
        str6='#binning '
        str7='#filter '
        str8='#interval '
        str9='#quitat '
        str10='#chain '
        s=''
        s2=''
        seq_ra=(str(int(ra1)),'h',str(int(ra2)),'m',str(ra3),'s')
        seq_dec=(str(int(dec1)),'d',str(abs(int(dec2))),'m',str(abs(dec3)),'s')
        s.join( seq_ra )
        s2.join( seq_dec )
        #print(s.join( seq_ra ),s2.join( seq_dec ))
        c = SkyCoord(s.join( seq_ra ),s2.join( seq_dec ),frame='icrs')
        #print(c.dec.degree,c.dec.radian)
        for i in range(1,2):
            out.write(str00 + '\n')
        out.write(str00 + ' ' + str(name).replace('_2','') + '\n')
        out.write(str00 + '\n')
        out.write(str00 + ' ' + str(gaia_id_target) + '\n')
        for i in range(1,2):
            out.write(str00 + '\n')
        out.write(str1 + str(hour0) + ':' + str(minute0) + '\n')
        out.write(str22)
        out.write('#chill -60\n')
        out.write(str2 + str(waitlimit) + '\n')
        out.write(str3)
        out.write(str33 + str(afinterval) +'\n')
        if autofocus is None:
            out.write(str00 + str4)
        out.write(str5 + str(count) + '\n')
        out.write(str6 + str(binning) + '\n')
        out.write(str7 + str(filt) + '\n')
        out.write(str8 + str(exptime) +'\n')
        out.write(';#TAG Donuts=on'+'\n')
        if int(dec1)==-0:
            #print('ici')
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n') #8*np.cos(c.dec.radian)
        if int(dec1)<0 and int(dec1)!=-0:
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(abs(dec1)))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        if int(dec1)>0:
            #print('la')
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '+' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        out.write(str00 + '\n')
        out.write(str9 + str(hour1) + ':' + str(minute1) + '\n' )
        if name_2 is None:
            out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
        else:
            out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
        out.write(str00 + '\n')

def target_offset(t_now,name,date_start,date_end,waitlimit,afinterval,autofocus,count,filt,exptime,ra1,ra2,ra3,dec1,dec2,dec3,name_2,Path):
    df = pd.read_csv(target_list_path,delimiter = ' ',index_col = False)
    idx_target = None
    if 'NOI-105435' in name:
        gaia_id_target = 6915818294923863040
    if 'Trappist' in name:
        idx_target = np.where((df['Sp_ID'] == 'Sp2306-0502'))[0]
        gaia_id_target = int(df['Gaia_ID'][idx_target].values)
    else:
        idx_target = np.where((df['Sp_ID'] == name.replace('_2','')))[0]
        gaia_id_target = int(df['Gaia_ID'][idx_target].values)
    if not idx_target:
        df2 = pd.read_csv(target_list_path,delimiter = ',',index_col = False)
        idx_target = np.where((df2['Sp_ID'] == name))[0]
        gaia_id_target = int(df['Gaia_ID'][idx_target].values)
    date_start=np.datetime64(date_start)
    date_end=np.datetime64(date_end)
    binning=1
    name_file=name + '.txt'
    hour0=date_start.astype(object).hour
    hour0='{:02d}'.format(int(hour0))
    minute0=date_start.astype(object).minute
    minute0='{:02d}'.format(int(minute0))
    hour1=date_end.astype(object).hour
    hour1='{:02d}'.format(int(hour1))
    minute1=date_end.astype(object).minute
    minute1='{:02d}'.format(int(minute1))
    with open(os.path.join(Path,str(t_now),'Obj_' + name_file),'w') as out:
        str00=';'
        str1='#waituntil 1, '
        str2='#waitinlimits '
        str3='#nopreview \n'
        str33=';#afinterval '
        str4='#autofocus \n'
        str5='#count '
        str6='#binning '
        str7='#filter '
        str8='#interval '
        str9='#quitat '
        str10='#chain '
        s=''
        s2=''
        seq_ra=(str(int(ra1)),'h',str(int(ra2)),'m',str(ra3),'s')
        seq_dec=(str(int(dec1)),'d',str(abs(int(dec2))),'m',str(abs(dec3)),'s')
        s.join( seq_ra )
        s2.join( seq_dec )
        #print(s.join( seq_ra ),s2.join( seq_dec ))
        c = SkyCoord(s.join( seq_ra ),s2.join( seq_dec ),frame='icrs')
        #print(c.dec.degree,c.dec.radian)
        teldeg = 2
        skyarcmin = (teldeg/60)/(np.cos(c.dec.degree*np.pi/180))*60
        ra_deg = c.ra.degree + skyarcmin/60
        dec_deg = c.dec.degree + teldeg/60
        c = SkyCoord(ra = ra_deg *u.degree, dec = dec_deg*u.degree)
        ra1 = c.ra.hms.h
        ra2 = c.ra.hms.m
        ra3 = c.ra.hms.s
        dec1 = c.dec.dms.d
        dec2 = c.dec.dms.m
        dec3 = c.dec.dms.s
        #print('ICI TOI' ,c.ra.hms,c.dec.dms,ra1,ra2,ra3,dec1,dec2,dec3)
        for i in range(1,2):
            out.write(str00 + '\n')
        out.write(str00 + ' ' + str(name.replace('_2','')) + '\n')
        out.write(str00 + '\n')
        out.write(str00 + ' ' + str(gaia_id_target) + '\n')
        for i in range(1,2):
            out.write(str00 + '\n')
        out.write(str1 + str(hour0) + ':' + str(minute0) + '\n')
        out.write('#chill -60\n')
        out.write(str2 + str(waitlimit) + '\n')
        out.write(str3)
        out.write(str33 + str(afinterval) +'\n')
        if autofocus is None:
            out.write(str00 + str4)
        out.write(str5 + str(count) + '\n')
        out.write(str6 + str(binning) + '\n')
        out.write(str7 + str(filt) + '\n')
        out.write(str8 + str(exptime) +'\n')
        out.write('#TAG Donuts=on'+'\n')
        if int(dec1)==-0:
            #print('ici')
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n') #8*np.cos(c.dec.radian)
        if int(dec1)<0 and int(dec1)!=-0:
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(abs(dec1)))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        if int(dec1)>0:
            #print('la')
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '+' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        out.write(str00 + '\n')
        out.write(str9 + str(hour1) + ':' + str(minute1) + '\n' )
        if name_2 is None:
            out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
        else:
            out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
        out.write(str00 + '\n')

def first_target(t_now,name,date_start,date_end,waitlimit,afinterval,autofocus,count,filt,exptime,ra1,ra2,ra3,dec1,dec2,dec3,name_2,Path,telescope=None):
    df = pd.read_csv(target_list_path,delimiter = ' ',index_col = False)
    #print(name)
    idx_target = None
    if 'Trappist' in name:
        idx_target = np.where((df['Sp_ID'] == 'Sp2306-0502'))[0]
        gaia_id_target = int(df['Gaia_ID'][idx_target].values)
    elif 'Sp' in name:
        df2 = pd.read_csv(target_list_path,delimiter=' ',index_col=None)
        #idx_target = np.where((df['spc'] == name.replace('_2','')))[0]
        idx_target = np.where((df2['Sp_ID'] == name.replace('_2','')))[0]
        gaia_id_target = df2['Gaia_ID'][int(idx_target)] #int(df['gaia'][idx_target].values)
    if idx_target is None:
        df2 = pd.read_csv('target_list_special.txt',delimiter=' ',index_col=None)
        try:
            idx_target = np.where((df2['Sp_ID'] == name))[0]
            gaia_id_target = df2['Gaia_ID'][int(idx_target)]
        except TypeError:
            gaia_id_target = 'None'

    date_start=np.datetime64(date_start)
    date_end=np.datetime64(date_end)
    binning=1
    name_file=name + '.txt'
    hour0=date_start.astype(object).hour
    hour0='{:02d}'.format(int(hour0))
    minute0=date_start.astype(object).minute
    minute0='{:02d}'.format(int(minute0))
    hour1=date_end.astype(object).hour
    hour1='{:02d}'.format(int(hour1))
    minute1=date_end.astype(object).minute
    minute1='{:02d}'.format(int(minute1))
    with open(os.path.join(Path,str(t_now),'Obj_' + name_file),'w') as out:
        str00=';'
        str1='#waituntil 1, '
        str2='#waitinlimits '
        str22='#domeopen\n'
        str3='#nopreview \n'
        str33=';#afinterval '
        str4='#autofocus \n'
        str5='#count '
        str6='#binning '
        str7='#filter '
        str8='#interval '
        str9='#quitat '
        str10='#chain '
        s=''
        s2=''
        seq_ra=(str(int(ra1)),'h',str(int(ra2)),'m',str(ra3),'s')
        seq_dec=(str(int(dec1)),'d',str(abs(int(dec2))),'m',str(abs(dec3)),'s')
        s.join( seq_ra )
        s2.join( seq_dec )
        #print(s.join( seq_ra ),s2.join( seq_dec ))
        c = SkyCoord(s.join( seq_ra ),s2.join( seq_dec ),frame='icrs')
        #print(c.dec.degree,c.dec.radian)
        for i in range(1,2):
            out.write(str00 + '\n')
        out.write(str00 + ' ' + str(name.replace('_2','')) + '\n')
        out.write(str00 + '\n')
        out.write(str00 + ' ' + str(gaia_id_target) + '\n')
        for i in range(1,2):
            out.write(str00 + '\n')
        out.write(str1 + pd.to_datetime(str(date_start)).strftime('%Y/%m/%d %H:%M:%S') + '\n')#+ str(hour0) + ':' + str(minute0) + '\n')
        out.write(str22)
        if telescope == 'Artemis':
            out.write('#chill -55\n')
        else:
            out.write('#chill -60\n')
        out.write(str2 + str(waitlimit) + '\n')
        out.write(str3)
        out.write(str33 + str(afinterval) +'\n')
        if autofocus is None:
            out.write(str00 + str4)
        out.write(str5 + str(count) + '\n')
        out.write(str6 + str(binning) + '\n')
        out.write(str7 + str(filt) + '\n')
        out.write(str8 + str(exptime) +'\n')
        out.write('#TAG Donuts=on'+'\n')
        if int(dec1)==-0:
            #print('ici')
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n') #8*np.cos(c.dec.radian)
        if int(dec1)<0 and int(dec1)!=-0:
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(abs(dec1)))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        if int(dec1)>0:
            #print('la')
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '+' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        out.write(str00 + '\n')
        out.write(str9 + str(hour1) + ':' + str(minute1) + '\n' )
        if name_2 is None:
            out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
        else:
            out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
        out.write(str00 + '\n')

def first_target_offset(t_now,name,date_start,date_end,waitlimit,afinterval,autofocus,count,filt,exptime,ra1,ra2,ra3,dec1,dec2,dec3,name_2,Path):
    date_start=np.datetime64(date_start)
    date_end=np.datetime64(date_end)
    binning=1
    name_file=name + '.txt'
    hour0=date_start.astype(object).hour
    hour0='{:02d}'.format(int(hour0))
    minute0=date_start.astype(object).minute
    minute0='{:02d}'.format(int(minute0))
    hour1=date_end.astype(object).hour
    hour1='{:02d}'.format(int(hour1))
    minute1=date_end.astype(object).minute
    minute1='{:02d}'.format(int(minute1))
    with open(os.path.join(Path,str(t_now),'Obj_' + name_file),'w') as out:
        str00=';'
        str1='#waituntil 1, '
        str2='#waitinlimits '
        str22='#domeopen\n'
        str3='#nopreview \n'
        str33=';#afinterval '
        str4='#autofocus \n'
        str5='#count '
        str6='#binning '
        str7='#filter '
        str8='#interval '
        str9='#quitat '
        str10='#chain '
        s=''
        s2=''
        seq_ra=(str(int(ra1)),'h',str(int(ra2)),'m',str(ra3),'s')
        seq_dec=(str(int(dec1)),'d',str(abs(int(dec2))),'m',str(abs(dec3)),'s')
        s.join( seq_ra )
        s2.join( seq_dec )
        #print(s.join( seq_ra ),s2.join( seq_dec ))
        c = SkyCoord(s.join( seq_ra ),s2.join( seq_dec ),frame='icrs')
        #print(c.dec.degree,c.dec.radian)
        teldeg = 2
        skyarcmin = (teldeg/60)/(np.cos(c.dec.degree*np.pi/180))*60
        ra_deg = c.ra.degree + skyarcmin/60
        dec_deg = c.dec.degree + teldeg/60
        c = SkyCoord(ra = ra_deg *u.degree, dec = dec_deg*u.degree)
        ra1 = c.ra.hms.h
        ra2 = c.ra.hms.m
        ra3 = c.ra.hms.s
        dec1 = c.dec.dms.d
        dec2 = c.dec.dms.m
        dec3 = c.dec.dms.s
        #print('ICI TOI first target' ,c.ra.hms,c.dec.dms,ra1,ra2,ra3,dec1,dec2,dec3)
        for i in range(1,2):
            out.write(str00 + '\n')
        out.write(str00 + ' ' + str(name) + '\n')
        for i in range(1,2):
            out.write(str00 + '\n')
        out.write(str1 + pd.to_datetime(str(date_start)).strftime('%Y/%m/%d %H:%M:%S') + '\n')#+ str(hour0) + ':' + str(minute0) + '\n')
        out.write(str22)
        out.write('#chill -60\n')
        out.write(str2 + str(waitlimit) + '\n')
        out.write(str3)
        out.write(str33 + str(afinterval) +'\n')
        if autofocus is None:
            out.write(str00 + str4)
        out.write(str5 + str(count) + '\n')
        out.write(str6 + str(binning) + '\n')
        out.write(str7 + str(filt) + '\n')
        out.write(str8 + str(exptime) +'\n')
        out.write('#TAG Donuts=on'+'\n')
        if int(dec1)==-0:
            #print('ici')
            out.write(name + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n') #8*np.cos(c.dec.radian)
        if int(dec1)<0 and int(dec1)!=-0:
            out.write(name + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(abs(dec1)))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        if int(dec1)>0:
            #print('la')
            out.write(name + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '+' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        out.write(str00 + '\n')
        out.write(str9 + str(hour1) + ':' + str(minute1) + '\n' )
        if name_2 is None:
            out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
        else:
            out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
        out.write(str00 + '\n')

def first_target_no_DONUTS(t_now,name,date_start,date_end,waitlimit,afinterval,autofocus,count,filt,exptime,ra1,ra2,ra3,dec1,dec2,dec3,name_2,Path):
    date_start=np.datetime64(date_start)
    date_end=np.datetime64(date_end)
    binning=1
    name_file=name + '.txt'
    hour0=date_start.astype(object).hour
    hour0='{:02d}'.format(int(hour0))
    minute0=date_start.astype(object).minute
    minute0='{:02d}'.format(int(minute0))
    hour1=date_end.astype(object).hour
    hour1='{:02d}'.format(int(hour1))
    minute1=date_end.astype(object).minute
    minute1='{:02d}'.format(int(minute1))
    with open(os.path.join(Path,str(t_now),'Obj_' + name_file),'w') as out:
        str00=';'
        str1='#waituntil 1, '
        str2='#waitinlimits '
        str22='#domeopen\n'
        str3='#nopreview \n'
        str33=';#afinterval '
        str4='#autofocus \n'
        str5='#count '
        str6='#binning '
        str7='#filter '
        str8='#interval '
        str9='#quitat '
        str10='#chain '
        s=''
        s2=''
        seq_ra=(str(int(ra1)),'h',str(int(ra2)),'m',str(ra3),'s')
        seq_dec=(str(int(dec1)),'d',str(abs(int(dec2))),'m',str(abs(dec3)),'s')
        s.join( seq_ra )
        s2.join( seq_dec )
        #print(s.join( seq_ra ),s2.join( seq_dec ))
        c = SkyCoord(s.join( seq_ra ),s2.join( seq_dec ),frame='icrs')
        #print(c.dec.degree,c.dec.radian)
        for i in range(1,2):
            out.write(str00 + '\n')
        out.write(str00 + ' ' + str(name) + '\n')
        for i in range(1,2):
            out.write(str00 + '\n')
        out.write(str1 + pd.to_datetime(str(date_start)).strftime('%Y/%m/%d %H:%M:%S') + '\n')#+ str(hour0) + ':' + str(minute0) + '\n')
        out.write(str22)
        out.write('#chill -60\n')
        out.write(str2 + str(waitlimit) + '\n')
        out.write(str3)
        out.write(str33 + str(afinterval) +'\n')
        if autofocus is None:
            out.write(str00 + str4)
        out.write(str5 + str(count) + '\n')
        out.write(str6 + str(binning) + '\n')
        out.write(str7 + str(filt) + '\n')
        out.write(str8 + str(exptime) +'\n')
        out.write(';#TAG Donuts=on'+'\n')
        if int(dec1)==-0:
            #print('ici')
            out.write(name + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n') #8*np.cos(c.dec.radian)
        if int(dec1)<0 and int(dec1)!=-0:
            out.write(name + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(abs(dec1)))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        if int(dec1)>0:
            #print('la')
            out.write(name + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '+' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        out.write(str00 + '\n')
        out.write(str9 + str(hour1) + ':' + str(minute1) + '\n' )
        if name_2 is None:
            out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
        else:
            out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
        out.write(str00 + '\n')

def flatdawn(t_now,date_end,sun_rise,Path):
    date_end=np.datetime64(date_end)
    hour0=date_end.astype(object).hour
    hour0='{:02d}'.format(int(hour0))
    minute0=date_end.astype(object).minute
    minute0='{:02d}'.format(int(minute0))
    sun_rise=np.datetime64(sun_rise)
    hour1=sun_rise.astype(object).hour
    hour1='{:02d}'.format(int(hour1))
    minute1=sun_rise.astype(object).minute
    minute1='{:02d}'.format(int(minute1))
    with open(os.path.join(Path,str(t_now),'Cal_flatdawn.txt'),'w') as out:
        str00=';'
        str1='#waituntil 1, '
        str3='#dawnflats '
        str4='Cal_flatexo.txt \n'
        str44='#quitat '
        str5='#chain '
        str6='Cal_biasdark.txt \n'
        out.write(str1 + str(hour0) + ':' + str(minute0) + '\n')
        out.write(str3 + str4)
        out.write(str44 + str(hour1) + ':' + str(minute1) + '\n')
        out.write(str5 + str6)

def flatdawn_artemis(t_now,date_end,sun_rise,Path):
    date_end=np.datetime64(date_end)
    hour0=date_end.astype(object).hour
    hour0='{:02d}'.format(int(hour0))
    minute0=date_end.astype(object).minute
    minute0='{:02d}'.format(int(minute0))
    sun_rise=np.datetime64(sun_rise)
    hour1=sun_rise.astype(object).hour
    hour1='{:02d}'.format(int(hour1))
    minute1=sun_rise.astype(object).minute
    minute1='{:02d}'.format(int(minute1))
    with open(os.path.join(Path,str(t_now),'Cal_flatdawn.txt'),'w') as out:
        str00=';'
        str1='#waituntil 1, '
        str3='#dawnflats '
        str4='Cal_flatexo_morning.txt \n'
        str44='#quitat '
        str5='#chain '
        str6='Cal_biasdark.txt \n'
        out.write(str1 + str(hour0) + ':' + str(minute0) + '\n')
        out.write(str3 + str4)
        out.write(str44 + str(hour1) + ':' + str(minute1) + '\n')
        out.write(str5 + str6)

def flatdawn_no_flats(t_now,date_end,sun_rise,Path):
    date_end=np.datetime64(date_end)
    hour0=date_end.astype(object).hour
    hour0='{:02d}'.format(int(hour0))
    minute0=date_end.astype(object).minute
    minute0='{:02d}'.format(int(minute0))
    sun_rise=np.datetime64(sun_rise)
    hour1=sun_rise.astype(object).hour
    hour1='{:02d}'.format(int(hour1))
    minute1=sun_rise.astype(object).minute
    minute1='{:02d}'.format(int(minute1))
    with open(os.path.join(Path,str(t_now),'Cal_flatdawn.txt'),'w') as out:
        str00=';'
        str1='#waituntil 1, '
        str3=';#dawnflats '
        str4=';Cal_flatexo.txt \n'
        str44=';#quitat '
        str5='#chain '
        str6='Cal_biasdark.txt \n'
        out.write(str1 + str(hour0) + ':' + str(minute0) + '\n')
        out.write(str3 + str4)
        out.write(str44 + str(hour1) + ':' + str(minute1) + '\n')
        out.write(str5 + str6)

def flatexo_gany(Path,t_now,filt,nbOIII=None,nbHa=None,nbSII=None,nbz=None,nbr=None,nbi=None,nbg=None,nbIz=None,nbExo=None,nbClear=None):
    str00=';'
    if nbOIII is None:
        nbOIII=3
    if nbHa is None:
        nbHa=3
    if nbSII is None:
        nbSII=3
    if nbz is None:
        nbz=3
    if nbr is None:
        nbr=3
    if nbi is None:
        nbi=3
    if nbg is None:
        nbg=3
    if nbIz is None:
        nbIz=3
    if nbExo is None:
        nbExo=3
    if nbClear is None:
        nbClear=3
    with open(os.path.join(Path,str(t_now),'Cal_flatexo.txt'),'w') as out:
        if filt.find('OIII',1)!=-1:
            out.write(str(nbOIII) + ',' + 'OIII' + ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbOIII) + ',' + 'OIII' + ',' + '1' + '\n')
        if filt.find('Ha',1)!=-1:
            out.write(str(nbHa) + ',' + 'Ha' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbHa) + ',' + 'Ha' ',' + '1' + '\n')
        if filt.find('SII',1)!=-1:
            out.write(str(nbSII) + ',' + 'SII' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbSII) + ',' + 'SII' + ',' + '1' + '\n')
        if filt.find('z\'',1)!=-1 and filt.find('I+z',1)==-1:
            out.write(str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        if filt.find('r\'',1)!=-1:
            out.write(str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        if filt.find('i\'',1)!=-1:
            out.write(str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        if filt.find('g\'',1)!=-1:
            out.write(str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        if filt.find('I+z',1)!=-1:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        if filt.find('Exo',1)!=-1:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        if filt.find('Clear',1)!=-1:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')

    # out.write(';5,Clear,1' + '\n')
    # out.write(';5,Clear1,1' + '\n')
    # out.write(';5,Clear2,1' + '\n')

def flatexo_calli(Path,t_now,filt,nbu=None,nbB=None,nbz=None,nbV=None,nbr=None,nbi=None,nbg=None,nbIz=None,nbExo=None,nbClear=None):#u=None,nbu=None,nbr=None,nbz=None,nbg=None,nbi=None,nbIz=None,nbExo=None):
    str00=';'
    if nbu is None:
        nbu=3
    if nbB is None:
        nbB=3
    if nbz is None:
        nbz=3
    if nbV is None:
        nbV=3
    if nbr is None:
        nbr=3
    if nbi is None:
        nbi=3
    if nbg is None:
        nbg=3
    if nbIz is None:
        nbIz=3
    if nbExo is None:
        nbExo=3
    if nbClear is None:
        nbClear=3
    with open(os.path.join(Path,str(t_now),'Cal_flatexo.txt'),'w') as out:
        if filt.find('u\'',1)!=-1:
            out.write(str(nbu) + ',' + 'u\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbu) + ',' + 'u\'' ',' + '1' + '\n')
        if filt.find('B',1)!=-1:
            out.write(str(nbB) + ',' + 'B' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbB) + ',' + 'B' + ',' + '1' + '\n')
        if filt.find('z\'',1)!=-1 and filt.find('I+z',1)==-1:
            out.write(str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        if filt.find('V',1)!=-1:
            out.write(str(nbV) + ',' + 'V' + ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbV) + ',' + 'V' + ',' + '1' + '\n')
        if filt.find('r\'',1)!=-1:
            out.write(str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        if filt.find('i\'',1)!=-1:
            out.write(str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        if filt.find('g\'',1)!=-1:
            out.write(str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        if filt.find('I+z',1)!=-1:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        if filt.find('Exo',1)!=-1:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        if filt.find('Clear',1)!=-1:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')

def flatexo_euro(Path,t_now,filt,nbRc=None,nbB=None,nbz=None,nbV=None,nbr=None,nbi=None,nbg=None,nbIz=None,nbExo=None,nbClear=None):#u=None,nbu=None,nbr=None,nbz=None,nbg=None,nbi=None,nbIz=None,nbExo=None):
    str00=';'
    if nbRc is None:
        nbRc=3
    if nbB is None:
        nbB=3
    if nbz is None:
        nbz=3
    if nbV is None:
        nbV=3
    if nbr is None:
        nbr=3
    if nbi is None:
        nbi=3
    if nbg is None:
        nbg=3
    if nbIz is None:
        nbIz=3
    if nbExo is None:
        nbExo=3
    if nbClear is None:
        nbClear=3
    with open(os.path.join(Path,str(t_now),'Cal_flatexo.txt'),'w') as out:
        if filt.find('Rc',1)!=-1:
            out.write(str(nbRc) + ',' + 'Rc' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbRc) + ',' + 'Rc' ',' + '1' + '\n')
        if filt.find('B',1)!=-1:
            out.write(str(nbB) + ',' + 'B' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbB) + ',' + 'B' + ',' + '1' + '\n')
        if filt.find('z\'',1)!=-1 and filt.find('I+z',1)==-1:
            out.write(str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        if filt.find('V',1)!=-1:
            out.write(str(nbV) + ',' + 'V' + ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbV) + ',' + 'V' + ',' + '1' + '\n')
        if filt.find('r\'',1)!=-1:
            out.write(str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        if filt.find('i\'',1)!=-1:
            out.write(str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        if filt.find('g\'',1)!=-1:
            out.write(str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        if filt.find('I+z',1)!=-1:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        if filt.find('Exo',1)!=-1:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        if filt.find('Clear',1)!=-1:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')

def flatexo_saintex(Path,t_now,filt,nbu=None,nbz=None,nbr=None,nbi=None,nbg=None,nbIz=None,nbExo=None,nbClear=None):#u=None,nbu=None,nbr=None,nbz=None,nbg=None,nbi=None,nbIz=None,nbExo=None):
    str00=';'
    if nbu is None:
        nbu=3
    if nbz is None:
        nbz=3
    if nbr is None:
        nbr=3
    if nbi is None:
        nbi=3
    if nbg is None:
        nbg=3
    if nbIz is None:
        nbIz=3
    if nbExo is None:
        nbExo=3
    if nbClear is None:
        nbClear=3
    with open(os.path.join(Path,str(t_now),'Cal_flatexo.txt'),'w') as out:
        if filt.find('u\'',1)!=-1:
            out.write(str(nbu) + ',' + 'u\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbu) + ',' + 'u\'' ',' + '1' + '\n')
        if filt.find('z\'',1)!=-1 and filt.find('I+z',1)==-1:
            out.write(str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        if filt.find('r\'',1)!=-1:
            out.write(str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        if filt.find('i\'',1)!=-1:
            out.write(str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        if filt.find('g\'',1)!=-1:
            out.write(str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        if filt.find('I+z',1)!=-1:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        if filt.find('Exo',1)!=-1:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        if filt.find('Clear',1)!=-1:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')

def flatexo_io(Path,t_now,filt,nbu=None,nbHa=None,nbRc=None,nbz=None,nbr=None,nbi=None,nbg=None,nbIz=None,nbExo=None,nbClear=None):
    str00=';'
    if nbu is None:
        nbu=3
    if nbHa is None:
        nbHa=3
    if nbRc is None:
        nbu=Rc
    if nbz is None:
        nbz=3
    if nbr is None:
        nbr=3
    if nbi is None:
        nbi=3
    if nbg is None:
        nbg=3
    if nbIz is None:
        nbIz=3
    if nbExo is None:
        nbExo=3
    if nbClear is None:
        nbClear=3
    with open(os.path.join(Path,str(t_now),'Cal_flatexo.txt'),'w') as out:
        if filt.find('u\'',1)!=-1:
            out.write(str(nbu) + ',' + 'u\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbu) + ',' + 'u\'' ',' + '1' + '\n')
        if filt.find('Ha',1)!=-1:
            out.write(str(nbHa) + ',' + 'Ha' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbHa) + ',' + 'Ha' + ',' + '1' + '\n')
        if filt.find('Rc',1)!=-1:
            out.write(str(nbRc) + ',' + 'Rc' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbRc) + ',' + 'Rc' + ',' + '1' + '\n')
        if filt.find('z\'',1)!=-1 and filt.find('I+z',1)==-1:
            out.write(str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        if filt.find('r\'',1)!=-1:
            out.write(str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        if filt.find('i\'',1)!=-1:
            out.write(str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        if filt.find('g\'',1)!=-1:
            out.write(str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        if filt.find('I+z',1)!=-1:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        if filt.find('Exo',1)!=-1:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        if filt.find('Clear',1)!=-1:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')

def flatexo_artemis_morning(Path,t_now,filt,nbu=None,nbz=None,nbr=None,nbi=None,nbg=None,nbIz=None,nbExo=None,nbClear=None):#u=None,nbu=None,nbr=None,nbz=None,nbg=None,nbi=None,nbIz=None,nbExo=None):
    str00=';'
    if nbu is None:
        nbu=7
    if nbz is None:
        nbz=7
    if nbr is None:
        nbr=7
    if nbi is None:
        nbi=7
    if nbg is None:
        nbg=7
    if nbIz is None:
        nbIz=7
    if nbExo is None:
        nbExo=7
    if nbClear is None:
        nbClear=7
    with open(os.path.join(Path,str(t_now),'Cal_flatexo_morning.txt'),'w') as out:
        if filt.find('u',1)!=-1:
            out.write(str(nbu) + ',' + 'u' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbu) + ',' + 'u' ',' + '1' + '\n')
        if filt.find('z',1)!=-1 and filt.find('I+z',1)==-1:
            out.write(str(nbz) + ',' + 'z' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z' ',' + '1' + '\n')
        if filt.find('r',1)!=-1:
            out.write(str(nbr) + ',' + 'r' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r' ',' + '1' + '\n')
        if filt.find('i',1)!=-1:
            out.write(str(nbi) + ',' + 'i' + ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbi) + ',' + 'i' ',' + '1' + '\n')
        if filt.find('g',1)!=-1:
            out.write(str(nbg) + ',' + 'g' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g' + ',' + '1' + '\n')
        if filt.find('I+z',1)!=-1:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        if filt.find('Exo',1)!=-1:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        if filt.find('Clear',1)!=-1:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')

def flatexo_artemis_evening(Path,t_now,filt,nbu=None,nbHa=None,nbRc=None,nbz=None,nbr=None,nbi=None,nbg=None,nbIz=None,nbExo=None,nbClear=None):#u=None,nbu=None,nbr=None,nbz=None,nbg=None,nbi=None,nbIz=None,nbExo=None):
    str00=';'
    if nbu is None:
        nbu=7
    if nbz is None:
        nbz=7
    if nbr is None:
        nbr=7
    if nbi is None:
        nbi=7
    if nbg is None:
        nbg=7
    if nbIz is None:
        nbIz=7
    if nbExo is None:
        nbExo=7
    if nbClear is None:
        nbClear=7
    with open(os.path.join(Path,str(t_now),'Cal_flatexo_evening.txt'),'w') as out:
        if filt.find('u',1)!=-1:
            out.write(str(nbu) + ',' + 'u' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbu) + ',' + 'u' ',' + '1' + '\n')
        if filt.find('z',1)!=-1 and filt.find('I+z',1)==-1:
            out.write(str(nbz) + ',' + 'z' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z' ',' + '1' + '\n')
        if filt.find('r',1)!=-1:
            out.write(str(nbr) + ',' + 'r' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r' ',' + '1' + '\n')
        if filt.find('i',1)!=-1:
            out.write(str(nbi) + ',' + 'i' + ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbi) + ',' + 'i' ',' + '1' + '\n')
        if filt.find('g',1)!=-1:
            out.write(str(nbg) + ',' + 'g' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g' + ',' + '1' + '\n')
        if filt.find('I+z',1)!=-1:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        if filt.find('Exo',1)!=-1:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        if filt.find('Clear',1)!=-1:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')

def biasdark(t_now,Path):
    with open(os.path.join(Path,str(t_now),'Cal_biasdark.txt'),'w') as out:
        str00=';'
        str1='#domeclose \n'
        str2='#nopreview \n'
        str3=' == bias dark exoplanet =='
        str4='#count '
        str5='#binning '
        str6='#interval '
        str7='#dark \n'
        str8='#shutdown \n'
        str9='END'
        out.write(str1)
        out.write(str2)
        out.write(str00 + '\n')
        out.write(str00 + str3 + '\n')
        out.write(str00 + '\n')
        out.write(str4 + '9,9,9,9,9' + '\n')
        out.write(str5 + '1,1,1,1,1' + '\n')
        out.write(str6 + '0,15,30,60,120' + '\n')
        out.write(str7)
        out.write(str8)
        out.write(str00 + '\n')
        out.write(str00 + str9+ '\n')
        out.write(str00 + '\n')

def biasdark_comete(t_now,Path):
    with open(os.path.join(Path,str(t_now),'Cal_biasdark.txt'),'w') as out:
        str00=';'
        str1='#domeclose \n'
        str2='#nopreview \n'
        str3=' == bias dark exoplanet =='
        str4='#count '
        str5='#binning '
        str6='#interval '
        str7='#dark \n'
        str8='#shutdown \n'
        str9='END'
        out.write(str1)
        out.write(str2)
        out.write(str00 + '\n')
        out.write(str00 + str3 + '\n')
        out.write(str00 + '\n')
        out.write(str4 + '9,9,9,9,9,9,9' + '\n')
        out.write(str5 + '1,1,1,1,2,2,2' + '\n')
        out.write(str6 + '120,15,30,60,0,120,240' + '\n')
        out.write(str7)
        out.write(str8)
        out.write(str00 + '\n')
        out.write(str00 + str9+ '\n')
        out.write(str00 + '\n')

def shutdown(t_now,Path):
    with open(os.path.join(Path,str(t_now),'Cal_shutdown.txt'),'w') as out:
        str00=';'
        str1='#domeclose \n'
        str2='#shutdown \n'
        out.write(str00 + '\n')
        out.write(str1)
        out.write(str2)
        out.write(str00 + '\n')