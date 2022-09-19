from astropy import units as u
import numpy as np
import os
import pandas as pd
from datetime import date, timedelta, datetime
from astroplan import Observer
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import gspread
from SPOCK import path_spock
from oauth2client.service_account import ServiceAccountCredentials
from SPOCK import path_credential_json, target_list_from_stargate_path

startup_time = []
hour = []
minute = []
scope = ['https://spreadsheets.google.com/feeds',
         'https://www.googleapis.com/auth/drive']
creds = ServiceAccountCredentials.from_json_keyfile_name(path_credential_json, scope)
client = gspread.authorize(creds)

sh = client.open('SPECULOOS WG6')
# Read stars lists
worksheet_follow_up = sh.worksheet("Annex_Targets_V1-PLANETS")
dataframe = pd.DataFrame(worksheet_follow_up.get_all_records())
target_table_spc_follow_up = dataframe.rename(columns={"sp_id": "Sp_ID", "gaia_dr2": "Gaia_ID",
                                                            "period": "P", "period_e": "P_err",
                                                            "duration": "W", "duration_e": "W_err",
                                                            "dec": "DEC", "ra": "RA",
                                                            "dec_err": "DEC_err", "ra_err": "RA_err"})
target_table_spc_follow_up['W'] /= 24
target_table_spc_follow_up['W_err'] /= 24
# Read follow up (planet candidates) list
worksheet_special = sh.worksheet("Annex_Targets_V2-STARS")
dataframe = pd.DataFrame(worksheet_special.get_all_records())
target_table_spc_special = dataframe.rename(columns={"spc": "Sp_ID", "gaia": "Gaia_ID", "dec": "DEC",
                                                  "ra": "RA", "dec_err": "DEC_err", "ra_err": "RA_err",
                                                  "mag_j": "J", "V_mag": "V"})


def charge_observatory(Name):
    observatories = []
    if 'SSO' in str(Name):
        location = EarthLocation.from_geodetic(-70.40300000000002*u.deg, -24.625199999999996*u.deg,2635.0000000009704*u.m)
        observatories.append(Observer(location=location, name="SSO", timezone="UTC"))

    if 'SNO' in str(Name):
        location_SNO = EarthLocation.from_geodetic(-16.50583131*u.deg, 28.2999988*u.deg, 2390*u.m)
        observatories.append(Observer(location=location_SNO, name="SNO", timezone="UTC"))

    if 'Saint-Ex' in str(Name):
        location_saintex = EarthLocation.from_geodetic(-115.48694444444445*u.deg, 31.029166666666665*u.deg, 2829.9999999997976*u.m)
        observatories.append(Observer(location=location_saintex, name="Saint-Ex", timezone="UTC"))

    if 'TS_La_Silla' in str(Name):
        location_TSlasilla = EarthLocation.from_geodetic(-70.73000000000002*u.deg, -29.25666666666666*u.deg, 2346.9999999988418*u.m)
        observatories.append(Observer(location=location_TSlasilla, name="TS_La_Silla", timezone="UTC"))

    if 'TN_Oukaimeden' in str(Name):
        location_TNOuka = EarthLocation.from_geodetic(-7.862263*u.deg,31.20516*u.deg, 2751*u.m)
        observatories.append(Observer(location=location_TNOuka, name="TN_Oukaimeden", timezone="UTC"))

    return observatories


def Path_txt_files(telescope):
    Path=os.path.join(path_spock + '/DATABASE',telescope,'Plans_by_date')
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


def startup(t_now,name,sun_set,date_start,Path,telescope):
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
    if telescope == 'Saint-Ex':
        with open(os.path.join(Path,str(t_now),'startup_' + datetime.strptime(t_now,'%Y-%m-%d').strftime('%m%d%Y') + '.txt'),'w') as out:
            str00=';\n'
            str0=' == Startup =='
            str1='#waituntil 1, '
            str2='#domeopen\n'
            str3='#chill '
            str33='-70\n'
            str4='#duskflats Cal_flatexo'+ '_' + datetime.strptime(t_now, '%Y-%m-%d').strftime('%m%d%Y') +'.txt\n'
            str5='#chain '
            str6='#quitat '
            observatory = charge_observatory('Saint-Ex')[0]
            start_time_saint_ex = Time((observatory.twilight_evening_civil(Time(t_now)+1,
                                            which='nearest').jd + Time(sun_set).jd)/2,format='jd').tt.datetime
            out.write(';' + str0 + '\n')
            for i in range(1,3):
                out.write(str00)
            out.write(str1 + pd.to_datetime(str(startup_time1)).strftime('%m/%d/%Y %H:%M:%S') + '\n')#+ str(hour00) + ':' + str(minute00) + '\n')
            out.write(str2)
            for i in range(1,2):
                out.write(str00)
            #out.write(str1 + str(hour2) + ':' + str(minute2) + '\n')
            out.write(str3 + str33)
            for i in range(1,2):
                out.write(str00)
            out.write(str1 + start_time_saint_ex.strftime('%m/%d/%Y %H:%M:%S')  + '\n')#str(hour1) + ':' + str(minute1) + '\n')
            out.write(str4)
            for i in range(1,2):
                out.write(str00)
            out.write(str6 + pd.to_datetime(str(date_start)).strftime('%m/%d/%Y %H:%M:%S') + '\n')#str(hour0) + ':' + str(minute0) + '\n')
            for i in range(1,2):
                out.write(str00)
            out.write(str5 + 'Obj_' + name + '_' + datetime.strptime(t_now, '%Y-%m-%d').strftime('%m%d%Y') + '.txt' + '\n')
            out.write(str00)

    else:
        if(telescope == 'Io') or (telescope == 'Europa')  or (telescope == 'Ganymede') or (telescope == 'Callisto'):
            with open(os.path.join(Path, str(t_now), 'default.txt'), 'w') as out:
                str00 = ';\n'
                str0 = ' == Startup =='
                str1 = '#waituntil 1, '
                str2 = '#domeopen\n'
                str3 = '#chill '
                str33 = '-60\n'
                str4 = '#duskflats Cal_flatexo.txt\n'
                str5 = '#chain '
                str6 = '#quitat '
                out.write(';' + str0 + '\n')
                for i in range(1, 3):
                    out.write(str00)
                out.write(str1 + pd.to_datetime(str(sun_set)).strftime(
                    '%Y/%m/%d %H:%M:%S') + '\n')  # + str(hour00) + ':' + str(minute00) + '\n')
                out.write(str2)
                for i in range(1, 2):
                    out.write(str00)
                # out.write(str1 + str(hour2) + ':' + str(minute2) + '\n')
                out.write(str3 + str33)
                for i in range(1, 2):
                    out.write(str00)
                out.write(str1 + pd.to_datetime(str(startup_time1)).strftime(
                    '%Y/%m/%d %H:%M:%S') + '\n')  # str(hour1) + ':' + str(minute1) + '\n')
                out.write(str4)
                for i in range(1, 2):
                    out.write(str00)
                out.write(str6 + pd.to_datetime(str(date_start)).strftime(
                    '%Y/%m/%d %H:%M:%S') + '\n')  # str(hour0) + ':' + str(minute0) + '\n')
                for i in range(1, 2):
                    out.write(str00)
                out.write(str5 + 'Obj_' + name + '.txt' + '\n')
                out.write(str00)
        else:
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
    with open(os.path.join(Path,str(t_now),'default.txt'),'w') as out:
        str00=';\n'
        str0=' == Startup =='
        str1='#waituntil 1, '
        str2='#domeopen\n'
        str3='#chill '
        str33='-60\n'
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


def target(t_now,name,date_start,date_end,waitlimit,afinterval,autofocus,count,filt,exptime,
           ra1, ra2, ra3, dec1, dec2, dec3, name_2, Path, telescope):
    df = pd.read_csv(target_list_from_stargate_path, delimiter=',', index_col=False)
    idx_target = None
    if 'Trappist' in name:
        idx_target = np.where((df['Sp_ID'] == 'Sp2306-0502'))[0]
        gaia_id_target = int(df['Gaia_ID'][idx_target].values)
    elif ('Sp' in name) and ('.0' not in name):
        if name != 'Sp1837+2030':
            df2 = pd.read_csv(target_list_from_stargate_path, delimiter=',', index_col=None)
            if name not in ["Sp1633-6808_2","Sp0933-4353_2","Sp1953+4424_2"]:
                if name.find("_zcut") != -1:
                    idx_target = np.where((df2['Sp_ID'] == name.replace('_zcut', '')))[0]
                else:
                    idx_target = np.where((df2['Sp_ID'] == name.replace('_2', '')))[0]
            else:
                idx_target = np.where((df2['Sp_ID'] == name))[0]
            gaia_id_target = df2['Gaia_ID'][int(idx_target)]
    if idx_target is None:
        df2 = pd.concat([target_table_spc_follow_up, target_table_spc_special], ignore_index=True)
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
    if telescope == 'Saint-Ex':
        with open(os.path.join(Path, str(t_now),'Obj_' + name + '_' +
                                                datetime.strptime(t_now, '%Y-%m-%d').strftime('%m%d%Y') + '.txt'),
                  'w') as out:
            str00 = ';'
            str1 = '#waituntil 1, '
            str2 = '#waitinlimits '
            str22 = '#domeopen\n'
            str3 = '#nopreview \n'
            str33 = ';#afinterval '
            str4 = '#autofocus \n'
            str5 = '#count '
            str6 = '#binning '
            str7 = '#filter '
            str8 = '#interval '
            str9 = '#quitat '
            str10 = '#chain '
            s = ''
            s2 = ''
            seq_ra = (str(int(ra1)), 'h', str(int(ra2)), 'm', str(ra3), 's')
            seq_dec = (str(int(dec1)), 'd', str(abs(int(dec2))), 'm', str(abs(dec3)), 's')
            s.join(seq_ra)
            s2.join(seq_dec)
            # print(s.join( seq_ra ),s2.join( seq_dec ))
            c = SkyCoord(s.join(seq_ra), s2.join(seq_dec), frame='icrs')
            # print(c.dec.degree,c.dec.radian)
            for i in range(1, 2):
                out.write(str00 + '\n')
            out.write(str00 + ' ' + str(name).replace('_2', '') + '\n')
            out.write(str00 + '\n')
            out.write(str00 + ' ' + str(gaia_id_target) + '\n')
            for i in range(1, 2):
                out.write(str00 + '\n')
            out.write(str1 + pd.to_datetime(str(date_start)).strftime('%m/%d/%Y %H:%M:%S') + '\n')
            out.write(str22)
            if telescope == 'Artemis':
                out.write('#chill -60\n')
            if telescope == 'Saint-Ex':
                out.write(';#chill -70\n')
            else:
                out.write('#chill -60\n')

            out.write(str2 + str(waitlimit) + '\n')
            out.write(str3)
            out.write(str33 + str(afinterval) + '\n')
            if not autofocus:
                out.write(str00 + str4)
            else:
                out.write(str4)
            out.write(str5 + str(count) + '\n')
            out.write(str6 + str(binning) + '\n')
            out.write(str7 + str(filt) + '\n')
            out.write(str8 + str(exptime) + '\n')
            out.write('#TAG Donuts=on' + '\n')
            if str(dec1) == '-0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                          '-' + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
            if str(dec1) == '0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t'
                         + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
            if int(dec1) < 0 and int(dec1) != -0:
                out.write(name.replace('_2', '') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str(
                    '{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str(
                    '{:02d}'.format(int(abs(dec1)))) \
                          + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str(
                    '{:05.2f}'.format(abs(dec3))) + '\n')
            if int(dec1) > 0:
                out.write(name.replace('_2', '') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str(
                    '{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '+' + str(
                    '{:02d}'.format(int(dec1))) \
                          + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str(
                    '{:05.2f}'.format(abs(dec3))) + '\n')
            out.write(str00 + '\n')
            out.write(str9 + pd.to_datetime(str(date_end)).strftime('%m/%d/%Y %H:%M:%S') + '\n')
            if name_2 is None:
                out.write(str10 + 'Cal_flatdawn' + '_' +
                          datetime.strptime(t_now, '%Y-%m-%d').strftime('%m%d%Y') + '.txt' + '\n')
            else:
                out.write(str10 + 'Obj_' + name_2 + '_' +
                          datetime.strptime(t_now, '%Y-%m-%d').strftime('%m%d%Y') + '.txt' + '\n')
            out.write(str00 + '\n')
    if telescope != 'Saint-Ex':
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
            s = ''
            s2 = ''
            seq_ra = (str(int(ra1)), 'h', str(int(ra2)), 'm', str(ra3), 's')
            seq_dec = (str(int(dec1)), 'd', str(abs(int(dec2))), 'm', str(abs(dec3)), 's')
            s.join(seq_ra)
            s2.join(seq_dec)
            c = SkyCoord(s.join(seq_ra), s2.join(seq_dec), frame='icrs')
            for i in range(1,2):
                out.write(str00 + '\n')
            out.write(str00 + ' ' + str(name).replace('_2', '') + '\n')
            out.write(str00 + '\n')
            out.write(str00 + ' ' + str(gaia_id_target) + '\n')
            for i in range(1, 2):
                out.write(str00 + '\n')
            out.write(str1 + str(hour0) + ':' + str(minute0) + '\n')
            if name == "dome_rot":
                out.write(';#domeopen\n')
            else:
                out.write(str22)
            if telescope == 'Artemis':
                out.write('#chill -60\n')
            else:
                if name == 'dome_rot':
                    out.write(';#chill -60\n')
                else:
                    out.write('#chill -60\n')
            if name == 'dome_rot':
                if telescope == "Callisto":
                    out.write(r"#dir C:\Users\SPIRIT\Documents\ACP Astronomy\Images\Dome_rot" + '\n')
                else:
                    out.write(r"#dir C:\Users\speculoos\Documents\ACP Astronomy\Images\Dome_rot" + '\n')
            if (name.find('Ch') != -1) or (name.find('ch') != -1) and gaia_id_target == 'None':
                d = Time(t_now).tt.datetime
                out.write(r"#dir C:\Users\speculoos\Documents\ACP Astronomy\Images\Chilean" + r"\ "[0] +
                          d.strftime("%Y%m%d") + '\n')
            out.write(str2 + str(waitlimit) + '\n')
            if name == 'dome_rot':
                out.write('#nopointing'+'\n')
            out.write(str3)
            out.write(str33 + str(afinterval) + '\n')
            if not autofocus:
                out.write(str00 + str4)
            else:
                out.write(str4)
            if name == 'dome_rot':
                out.write('#count 1' + '\n')
            else:
                out.write(str5 + str(count) + '\n')
            out.write(str6 + str(binning) + '\n')
            out.write(str7 + str(filt) + '\n')
            out.write(str8 + str(exptime) +'\n')
            if name == 'dome_rot':
                out.write('#TAG Donuts=off' + '\n')
            else:
                out.write('#TAG Donuts=on'+'\n')
            if str(dec1) == '-0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                          '-' + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
            if str(dec1) == '0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t'
                         + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
            if int(dec1) < 0 and str(dec1) != '-0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                          '-' + str('{:02d}'.format(int(abs(dec1)))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
            if int(dec1) > 0:
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                          '+' + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
            out.write(str00 + '\n')
            if name == 'dome_rot':
                out.write('#dir'+ '\n')
            out.write(str9 + str(hour1) + ':' + str(minute1) + '\n' )
            if name_2 is None:
                if telescope == 'Artemis':
                    out.write(str10 + 'Cal_flatdawn.txt' + '\n')
                else:
                    out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
            else:
                out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
            out.write(str00 + '\n')


def target_no_DONUTS(t_now,name,date_start,date_end,waitlimit,afinterval,autofocus,count,filt,exptime,
                     ra1,ra2,ra3,dec1,dec2,dec3,name_2,Path):
    df = pd.read_csv(target_list_from_stargate_path, delimiter=',', index_col=False)
    idx_target = None
    if 'Trappist' in name:
        idx_target = np.where((df['Sp_ID'] == 'Sp2306-0502'))[0]
        gaia_id_target = int(df['Gaia_ID'][idx_target].values)
    elif 'Sp' in name:
        df2 = pd.read_csv(target_list_from_stargate_path, delimiter=',', index_col=None)
        # idx_target = np.where((df['spc'] == name.replace('_2','')))[0]
        idx_target = np.where((df2['Sp_ID'] == name.replace('_2','')))[0]
        gaia_id_target = df2['Gaia_ID'][int(idx_target)]  # int(df['gaia'][idx_target].values)
    if idx_target is None:
        df2 = pd.concat([target_table_spc_follow_up, target_table_spc_special], ignore_index=True)
        idx_target = np.where((df2['Sp_ID'] == name))[0]
        gaia_id_target = df2['Gaia_ID'][int(idx_target)]

    date_start = np.datetime64(date_start)
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
        out.write(str33 + str(afinterval) + '\n')
        if not autofocus:
            out.write(str00 + str4)
        else:
            out.write(str4)
        out.write(str5 + str(count) + '\n')
        out.write(str6 + str(binning) + '\n')
        out.write(str7 + str(filt) + '\n')
        out.write(str8 + str(exptime) +'\n')
        out.write(';#TAG Donuts=on'+'\n')
        if str(dec1) == '-0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                          '-' + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
        if str(dec1) == '0.0':
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                      str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t'
                     + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
        if int(dec1) < 0 and int(dec1) != -0:
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                      str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                      '-' + str('{:02d}'.format(int(abs(dec1)))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        if int(dec1) > 0:
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                      str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                      '+' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        out.write(str00 + '\n')
        out.write(str9 + str(hour1) + ':' + str(minute1) + '\n' )
        if name_2 is None:
            out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
        else:
            out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
        out.write(str00 + '\n')


def target_offset(t_now,name,date_start,date_end,waitlimit,afinterval,autofocus,count,filt,exptime,
                  ra1,ra2,ra3,dec1,dec2,dec3,name_2,Path):
    df = pd.read_csv(target_list_from_stargate_path, delimiter=',', index_col=False)
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
        df2 = pd.concat([target_table_spc_follow_up, target_table_spc_special], ignore_index=True)
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
        out.write(str33 + str(afinterval) + '\n')
        if not autofocus:
            out.write(str00 + str4)
        else:
            out.write(str4)
        out.write(str5 + str(count) + '\n')
        out.write(str6 + str(binning) + '\n')
        out.write(str7 + str(filt) + '\n')
        out.write(str8 + str(exptime) +'\n')
        out.write('#TAG Donuts=on'+'\n')
        if str(dec1) == '-0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                          '-' + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
        if str(dec1) == '0.0':
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                      str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t'
                     + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
        if int(dec1) < 0 and int(dec1) != -0:
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                      str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                      '-' + str('{:02d}'.format(int(abs(dec1)))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        if int(dec1) > 0:
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                      str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                      '+' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        out.write(str00 + '\n')
        out.write(str9 + str(hour1) + ':' + str(minute1) + '\n' )
        if name_2 is None:
            out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
        else:
            out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
        out.write(str00 + '\n')


def first_target(t_now,name,date_start,date_end,waitlimit,afinterval,autofocus,count,filt,exptime,
                 ra1,ra2,ra3,dec1,dec2,dec3,name_2,Path,telescope):
    df = pd.read_csv(target_list_from_stargate_path, delimiter=',', index_col=False)
    idx_target = None
    if 'Trappist' in name:
        idx_target = np.where((df['Sp_ID'] == 'Sp2306-0502'))[0]
        gaia_id_target = int(df['Gaia_ID'][idx_target].values)
    elif 'Sp' in name:
        if name != 'Sp1837+2030':
            df2 = pd.read_csv(target_list_from_stargate_path, delimiter=',', index_col=None)
            if name not in ["Sp1633-6808_2","Sp0933-4353_2","Sp1953+4424_2"]:
                if name.find("_zcut") != -1:
                    idx_target = np.where((df2['Sp_ID'] == name.replace('_zcut', '')))[0]
                else:
                    idx_target = np.where((df2['Sp_ID'] == name.replace('_2', '')))[0]
            else:
                idx_target = np.where((df2['Sp_ID'] == name))[0]
            gaia_id_target = df2['Gaia_ID'][int(idx_target)]
    if idx_target is None:
        df2 = pd.concat([target_table_spc_special, target_table_spc_follow_up], ignore_index=True)
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
    if telescope == 'Saint-Ex':
        with open(os.path.join(Path, str(t_now),'Obj_' + name + '_' +
                                                datetime.strptime(t_now, '%Y-%m-%d').strftime('%m%d%Y') + '.txt'),
                  'w') as out:
            str00 = ';'
            str1 = '#waituntil 1, '
            str2 = '#waitinlimits '
            str22 = '#domeopen\n'
            str3 = '#nopreview \n'
            str33 = ';#afinterval '
            str4 = '#autofocus \n'
            str5 = '#count '
            str6 = '#binning '
            str7 = '#filter '
            str8 = '#interval '
            str9 = '#quitat '
            str10 = '#chain '
            s = ''
            s2 = ''
            seq_ra = (str(int(ra1)), 'h', str(int(ra2)), 'm', str(ra3), 's')
            seq_dec = (str(int(dec1)), 'd', str(abs(int(dec2))), 'm', str(abs(dec3)), 's')
            s.join(seq_ra)
            s2.join(seq_dec)
            # print(s.join( seq_ra ),s2.join( seq_dec ))
            c = SkyCoord(s.join(seq_ra), s2.join(seq_dec), frame='icrs')
            # print(c.dec.degree,c.dec.radian)
            for i in range(1, 2):
                out.write(str00 + '\n')
            out.write(str00 + ' ' + str(name.replace('_2', '')) + '\n')
            out.write(str00 + '\n')
            out.write(str00 + ' ' + str(gaia_id_target) + '\n')
            for i in range(1, 2):
                out.write(str00 + '\n')
            out.write(str1 + pd.to_datetime(str(date_start)).strftime('%m/%d/%Y %H:%M:%S') + '\n')  # + str(hour0) + ':' + str(minute0) + '\n')
            out.write(str22)
            if telescope == 'Artemis':
                out.write('#chill -60\n')
            if telescope == 'Saint-Ex':
                out.write(';#chill -70\n')
            else:
                out.write('#chill -60\n')
            out.write(str2 + str(waitlimit) + '\n')
            out.write(str3)
            out.write(str33 + str(afinterval) + '\n')
            if not autofocus:
                out.write(str00 + str4)
            else:
                out.write(str4)
            out.write(str5 + str(count) + '\n')
            out.write(str6 + str(binning) + '\n')
            out.write(str7 + str(filt) + '\n')
            out.write(str8 + str(exptime) + '\n')
            out.write('#TAG Donuts=on' + '\n')
            if str(dec1) == '-0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                          '-' + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
            if str(dec1) == '0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t'
                         + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
            if int(dec1) < 0 and int(dec1) != -0:
                out.write(name.replace('_2', '') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str(
                    '{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str(
                    '{:02d}'.format(int(abs(dec1)))) \
                          + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str(
                    '{:05.2f}'.format(abs(dec3))) + '\n')
            if int(dec1) > 0:
                out.write(name.replace('_2', '') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str(
                    '{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '+' + str(
                    '{:02d}'.format(int(dec1))) \
                          + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str(
                    '{:05.2f}'.format(abs(dec3))) + '\n')
            out.write(str00 + '\n')
            out.write(str9 + pd.to_datetime(str(date_end)).strftime('%m/%d/%Y %H:%M:%S') + '\n')
            if name_2 is None:
                out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
            else:
                out.write(str10 + 'Obj_' + name_2 + '_' +
                          datetime.strptime(t_now, '%Y-%m-%d').strftime('%m%d%Y') + '.txt' + '\n')
            out.write(str00 + '\n')

    if telescope != 'Saint-Ex':
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
            out.write(str1 + pd.to_datetime(str(date_start)).strftime('%Y/%m/%d %H:%M:%S') + '\n')  # + str(hour0) + ':' + str(minute0) + '\n')
            out.write(str22)
            if telescope == 'Artemis':
                out.write('#chill -60\n')
            if telescope == 'Saint-Ex':
                out.write(';#chill -70\n')
            else:
                out.write('#chill -60\n')
            if (name.find('Ch') != -1) or (name.find('ch') != -1) and gaia_id_target == 'None':
                d = Time(t_now).tt.datetime
                out.write(r"#dir C:\Users\speculoos\Documents\ACP Astronomy\Images\Chilean" + r"\ "[0] +
                          d.strftime("%Y%m%d") + '\n')
            out.write(str2 + str(waitlimit) + '\n')
            out.write(str3)
            out.write(str33 + str(afinterval) +'\n')
            if not autofocus:
                out.write(str00 + str4)
            else:
                out.write(str4)
            out.write(str5 + str(count) + '\n')
            out.write(str6 + str(binning) + '\n')
            out.write(str7 + str(filt) + '\n')
            out.write(str8 + str(exptime) +'\n')
            out.write('#TAG Donuts=on'+'\n')
            if str(dec1) == '-0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                          '-' + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
            if str(dec1) == '0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t'
                         + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
            if int(dec1)<0 and int(dec1)!=-0:
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                          '-' + str('{:02d}'.format(int(abs(dec1)))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
            if int(dec1)>0:
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                          '+' + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
            out.write(str00 + '\n')
            out.write(str9 + str(hour1) + ':' + str(minute1) + '\n' )
            if name_2 is None:
                if telescope == 'Artemis':
                    out.write(str10 + 'Cal_flatdawn.txt' + '\n')
                else:
                    out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
            else:
                out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
            out.write(str00 + '\n')


def first_target_offset(t_now,name,date_start,date_end,waitlimit,afinterval,autofocus,count,
                        filt,exptime,ra1,ra2,ra3,dec1,dec2,dec3,name_2,Path):
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
        if not autofocus:
            out.write(str00 + str4)
        else:
            out.write(str4)
        out.write(str5 + str(count) + '\n')
        out.write(str6 + str(binning) + '\n')
        out.write(str7 + str(filt) + '\n')
        out.write(str8 + str(exptime) +'\n')
        out.write('#TAG Donuts=on'+'\n')
        if str(dec1) == '-0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                          '-' + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
        if str(dec1) == '0.0':
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                      str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t'
                     + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
        if int(dec1)<0 and int(dec1)!=-0:
            out.write(name + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(abs(dec1)))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        if int(dec1)>0:
            out.write(name + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '+' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        out.write(str00 + '\n')
        out.write(str9 + str(hour1) + ':' + str(minute1) + '\n' )
        if name_2 is None:
            out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
        else:
            out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
        out.write(str00 + '\n')


def first_target_no_DONUTS(t_now,name,date_start,date_end,waitlimit,afinterval,autofocus,count,
                           filt,exptime,ra1,ra2,ra3,dec1,dec2,dec3,name_2,Path):
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
        if not autofocus:
            out.write(str00 + str4)
        else:
            out.write(str4)
        out.write(str5 + str(count) + '\n')
        out.write(str6 + str(binning) + '\n')
        out.write(str7 + str(filt) + '\n')
        out.write(str8 + str(exptime) +'\n')
        out.write(';#TAG Donuts=on'+'\n')
        if str(dec1) == '-0.0':
                out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                          str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' +
                          '-' + str('{:02d}'.format(int(dec1))) \
                 + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
        if str(dec1) == '0.0':
            out.write(name.replace('_2','') + '\t' + str('{:02d}'.format(int(ra1))) + ' ' +
                      str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t'
                     + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')  # 8*np.cos(c.dec.radian)
        if int(dec1) < 0 and int(dec1) != -0:
            out.write(name + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '-' + str('{:02d}'.format(int(abs(dec1)))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        if int(dec1) > 0:
            out.write(name + '\t' + str('{:02d}'.format(int(ra1))) + ' ' + str('{:02d}'.format(int(ra2))) + ' ' + str('{:05.2f}'.format(float(ra3))) + '\t' + '+' + str('{:02d}'.format(int(dec1))) \
             + ' ' + str('{:02d}'.format(int(abs(dec2)))) + ' ' + str('{:05.2f}'.format(abs(dec3))) + '\n')
        out.write(str00 + '\n')
        out.write(str9 + str(hour1) + ':' + str(minute1) + '\n' )
        if name_2 is None:
            out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
        else:
            out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
        out.write(str00 + '\n')


def flatdawn(t_now,date_end,sun_rise,Path,telescope):
    date_end_saint_ex = date_end
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
    if telescope == 'Saint-Ex':
        with open(os.path.join(Path, str(t_now), 'Cal_flatdawn' + '_' + datetime.strptime(t_now, '%Y-%m-%d').strftime('%m%d%Y') + '.txt'), 'w') as out:
            str00 = ';'
            str1 = '#waituntil 1, '
            str3 = '#dawnflats '
            str4 = 'Cal_flatexo'+ '_' + datetime.strptime(t_now, '%Y-%m-%d').strftime('%m%d%Y') +'.txt\n'
            str44 = '#quitat '
            str5 = '#chain '
            str6 = 'Cal_biasdark'+ '_' + datetime.strptime(t_now, '%Y-%m-%d').strftime('%m%d%Y') +'.txt\n'
            observatory = charge_observatory('Saint-Ex')[0]
            end_time_saint_ex = Time((observatory.twilight_morning_nautical(Time(t_now)+1,which='next').jd +
                                      observatory.twilight_morning_civil(Time(t_now)+1,which='next').jd)/2,
                                     format='jd').tt.datetime
            out.write(str1 + datetime.strptime(date_end_saint_ex, '%Y-%m-%d %H:%M:%S.%f').strftime('%m/%d/%Y %H:%M:%S') + '\n')
            out.write(str3 + str4)
            out.write(str44 + end_time_saint_ex.strftime('%m/%d/%Y %H:%M:%S')  + '\n')
            out.write(str5 + str6)
    else:
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


def flatexo_gany(Path,t_now,filt, nbOIII=None, nbHa=None, nbSII=None, nbz=None, nbr=None, nbi=None, nbg=None,
                 nbIz=None, nbExo=None, nbClear=None):
    str00=';'
    if nbOIII is None:
        nbOIII = 5
    if nbHa is None:
        nbHa = 5
    if nbSII is None:
        nbSII = 5
    if nbz is None:
        nbz = 5
    if nbr is None:
        nbr = 5
    if nbi is None:
        nbi = 5
    if nbg is None:
        nbg = 5
    if nbIz is None:
        nbIz=7
    if nbExo is None:
        nbExo = 5
    if nbClear is None:
        nbClear = 5
    with open(os.path.join(Path,str(t_now),'Cal_flatexo.txt'),'w') as out:
        if 'OIII' in filt:
            out.write(str(nbOIII) + ',' + 'OIII' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbOIII) + ',' + 'OIII' + ',' + '1' + '\n')

        if 'Ha' in filt:
            out.write(str(nbHa) + ',' + 'Ha' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbHa) + ',' + 'Ha' ',' + '1' + '\n')

        if 'SII' in filt:
            out.write(str(nbSII) + ',' + 'SII' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbSII) + ',' + 'SII' + ',' + '1' + '\n')

        if ('z' in filt) or ('z\'' in filt):
            out.write(str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z\'' ',' + '1' + '\n')

        if ('r' in filt) or ('r\'' in filt):
            out.write(str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r\'' ',' + '1' + '\n')

        if ('i' in filt) or ('i\'' in filt):
            out.write(str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbi) + ',' + 'i\'' ',' + '1' + '\n')

        if ('gg' in filt) or ('gg\'' in filt):
            out.write(str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g\'' ',' + '1' + '\n')

        if ('I+z' in filt) or ('I+z\'' in filt):
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')

        if 'Exo' in filt:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')

        if 'Clear' in filt:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')


def flatexo_calli(Path,t_now,filt, nbB=None, nbz=None, nbzcut=None, nbV=None, nbr=None, nbi=None, nbg=None,
                  nbIz=None, nbExo=None, nbClear=None, nbxYJ=None):  # u=None, nbu=None, nbr=None, nbz=None, nbg=None, nbi=None, nbIz=None, nbExo=None):
    str00=';'
    if nbB is None:
        nbB=5
    if nbz is None:
        nbz=5
    if nbzcut is None:
        nbzcut=5
    if nbV is None:
        nbV=5
    if nbr is None:
        nbr=5
    if nbi is None:
        nbi=5
    if nbg is None:
        nbg=5
    if nbIz is None:
        nbIz=7
    if nbExo is None:
        nbExo=5
    if nbClear is None:
        nbClear=20
    if nbxYJ is None:
        nbxYJ =20
    with open(os.path.join(Path,str(t_now),'Cal_flatexo.txt'),'w') as out:
        if ('B' in filt) or ('B\'' in filt):
            out.write(str(nbB) + ',' + 'B' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbB) + ',' + 'B' + ',' + '1' + '\n')

        if ('z' in filt) or ('z\'' in filt):
            out.write(str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z\'' ',' + '1' + '\n')

        if 'zcut' in filt:
            out.write(str(nbzcut) + ',' + 'zcut' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbzcut) + ',' + 'z\'' ',' + '1' + '\n')

        if ('V' in filt) or ('V\'' in filt):
            out.write(str(nbV) + ',' + 'V' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbV) + ',' + 'V' + ',' + '1' + '\n')

        if ('r' in filt) or ('r\'' in filt):
            out.write(str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r\'' ',' + '1' + '\n')

        if ('i' in filt) or ('i\'' in filt):
            out.write(str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbi) + ',' + 'i\'' ',' + '1' + '\n')

        if ('g' in filt) or ('g\'' in filt):
            out.write(str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g\'' ',' + '1' + '\n')

        if ('I+z' in filt) or ('I+z\'' in filt):
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')

        if 'Exo' in filt:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        if 'zYJ' in filt:
            out.write(str(nbxYJ) + ',' + 'zYJ' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbxYJ) + ',' + 'zYJ' + ',' + '1' + '\n')
        if 'Clear' in filt:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')



def flatexo_euro(Path,t_now,filt, nbRc=None, nbB=None, nbz=None, nbV=None, nbr=None,
                 nbi=None, nbg=None, nbIz=None, nbExo=None, nbClear=None):  # u=None, nbu=None, nbr=None, nbz=None, nbg=None, nbi=None, nbIz=None, nbExo=None):
    str00=';'
    if nbRc is None:
        nbRc = 5
    if nbB is None:
        nbB = 5
    if nbz is None:
        nbz = 5
    if nbV is None:
        nbV = 5
    if nbr is None:
        nbr = 5
    if nbi is None:
        nbi = 5
    if nbg is None:
        nbg = 5
    if nbIz is None:
        nbIz = 7
    if nbExo is None:
        nbExo = 5
    if nbClear is None:
        nbClear = 5
    with open(os.path.join(Path,str(t_now),'Cal_flatexo.txt'),'w') as out:
        if 'Rc' in filt:
            out.write(str(nbRc) + ',' + 'Rc' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbRc) + ',' + 'Rc' ',' + '1' + '\n')

        if ('BB' in filt) or ('BB\'' in filt):
            out.write(str(nbB) + ',' + 'B' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbB) + ',' + 'B' + ',' + '1' + '\n')

        if ('z' in filt) or ('z\'' in filt):
            out.write(str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z\'' ',' + '1' + '\n')

        if ('V' in filt) or ('V\'' in filt):
            out.write(str(nbV) + ',' + 'V' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbV) + ',' + 'V' + ',' + '1' + '\n')

        if ('r' in filt) or ('r\'' in filt):
            out.write(str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r\'' ',' + '1' + '\n')

        if ('i' in filt) or ('i\'' in filt):
            out.write(str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbi) + ',' + 'i\'' ',' + '1' + '\n')

        if ('g' in filt) or ('g\'' in filt):
            out.write(str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g\'' ',' + '1' + '\n')

        if ('I+z' in filt) or ('I+z\'' in filt):
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')

        if 'Exo' in filt:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')

        if 'Clear' in filt:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')


def flatexo_saintex(Path, t_now, filt, nbu=None, nbz=None, nbr=None, nbi=None, nbg=None,
                    nbIz=None, nbExo=None, nbClear=None):  # u=None, nbu=None, nbr=None, nbz=None, nbg=None, nbi=None, nbIz=None, nbExo=None):
    str00=';'
    if nbu is None:
        nbu=5
    if nbz is None:
        nbz=5
    if nbr is None:
        nbr=5
    if nbi is None:
        nbi=5
    if nbg is None:
        nbg=5
    if nbIz is None:
        nbIz=7
    if nbExo is None:
        nbExo=5
    if nbClear is None:
        nbClear=5
    with open(os.path.join(Path, str(t_now), 'Cal_flatexo' + '_' +
                                           datetime.strptime(t_now, '%Y-%m-%d').strftime('%m%d%Y') +
                                           '.txt'), 'w') as out:
        if ('u' in filt) or( 'u\'' in filt):
            out.write(str(nbu) + ',' + 'u\'' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbu) + ',' + 'u\'' ',' + '1' + '\n')

        if ('z' in filt) or('z\'' in filt):
            out.write(str(nbz) + ',' + 'z' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z\'' ',' + '1' + '\n')

        if ('r' in filt) or('r\'' in filt):
            out.write(str(nbr) + ',' + 'r\'' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r\'' ',' + '1' + '\n')

        if ('i' in filt) or('i\'' in filt):
            out.write(str(nbi) + ',' + 'i\'' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbi) + ',' + 'i\'' ',' + '1' + '\n')

        if ('g' in filt) or('g\'' in filt):
            out.write(str(nbg) + ',' + 'g\'' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g\'' + ',' + '1' + '\n')

        if ('I+z' in filt) or( 'I+z\'' in filt):
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')

        if 'Exo' in filt:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')

        if 'Clear' in filt:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')


def flatexo_io(Path, t_now, filt, nbu=None, nbHa=None, nbRc=None, nbz=None, nbr=None, nbi=None,
               nbg=None, nbIz=None, nbExo=None, nbClear=None, nbzcut=None):
    str00=';'
    if nbu is None:
        nbu=5
    if nbHa is None:
        nbHa=5
    if nbRc is None:
        nbRc=5
    if nbzcut is None:
        nbzcut=7
    if nbz is None:
        nbz=5
    if nbr is None:
        nbr=5
    if nbi is None:
        nbi=5
    if nbg is None:
        nbg=5
    if nbIz is None:
        nbIz=7
    if nbExo is None:
        nbExo=5
    if nbClear is None:
        nbClear=5
    with open(os.path.join(Path, str(t_now), 'Cal_flatexo.txt'), 'w') as out:
        if ('u' in filt) or ('u\'' in filt):
            out.write(str(nbu) + ',' + 'u' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbu) + ',' + 'u' ',' + '1' + '\n')

        if 'Ha' in filt:
            out.write(str(nbHa) + ',' + 'Ha' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbHa) + ',' + 'Ha' + ',' + '1' + '\n')

        if 'Rc' in filt:
            out.write(str(nbRc) + ',' + 'Rc' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbRc) + ',' + 'Rc' + ',' + '1' + '\n')
        if ('zcut' in filt):
            out.write(str(nbzcut) + ',' + 'zcut' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbzcut) + ',' + 'zcut' ',' + '1' + '\n')
        if ('z' in filt) or ('z\'' in filt):
            out.write(str(nbz) + ',' + 'z\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z\'' ',' + '1' + '\n')

        if ('r' in filt) or ('r\'' in filt):
            out.write(str(nbr) + ',' + 'r\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r\'' ',' + '1' + '\n')

        if ('i' in filt) or ('i\'' in filt):
            out.write(str(nbi) + ',' + 'i\'' ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbi) + ',' + 'i\'' ',' + '1' + '\n')

        if ('g' in filt) or ('g\'' in filt):
            out.write(str(nbg) + ',' + 'g\'' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g\'' ',' + '1' + '\n')

        if ('I+z' in filt) or ('I+z\'' in filt):
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')

        if 'Exo' in filt:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')

        if 'Clear' in filt:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')


def flatexo_artemis_morning(Path,t_now,filt, nbu=None, nbz=None, nbr=None, nbi=None, nbg=None,
                            nbIz=None, nbExo=None, nbClear=None, nbzcut=None):  # u=None, nbu=None, nbr=None, nbz=None, nbg=None, nbi=None, nbIz=None, nbExo=None):
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
    if nbzcut is None:
        nbzcut=7
    with open(os.path.join(Path, str(t_now), 'Cal_flatexo_morning.txt'), 'w') as out:
        if ('u' in filt) or ('u\'' in filt):
            out.write(str(nbu) + ',' + 'u' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbu) + ',' + 'u' ',' + '1' + '\n')

        if ('z' in filt) or ('z\'' in filt):
            out.write(str(nbz) + ',' + 'z' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z' ',' + '1' + '\n')

        if 'zcut' in filt:
            out.write(str(nbzcut) + ',' + 'zcut' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbzcut) + ',' + 'zcut' + ',' + '1' + '\n')

        if ('r' in filt) or ('r\'' in filt):
            out.write(str(nbr) + ',' + 'r' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r' ',' + '1' + '\n')

        if ('i' in filt) or ('i\'' in filt):
            out.write(str(nbi) + ',' + 'i' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbi) + ',' + 'i' ',' + '1' + '\n')

        if ('g' in filt) or ('g\'' in filt):
            out.write(str(nbg) + ',' + 'g' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g' + ',' + '1' + '\n')

        if ('I+z' in filt) or ('I+z\'' in filt):
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')

        if 'Exo' in filt:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')

        if 'Clear' in filt:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')




def flatexo_artemis_evening(Path,t_now,filt, nbu=None,
                            nbz=None, nbr=None, nbi=None, nbg=None, nbIz=None, nbExo=None, nbClear=None, nbzcut=None):  # u=None, nbu=None, nbr=None, nbz=None, nbg=None, nbi=None, nbIz=None, nbExo=None):
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
    if nbzcut is None:
        nbzcut=7
    with open(os.path.join(Path,str(t_now),'Cal_flatexo_evening.txt'),'w') as out:
        if ('u' in filt) or( 'u\'' in filt):
            out.write(str(nbu) + ',' + 'u' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbu) + ',' + 'u' ',' + '1' + '\n')

        if ('z' in filt) or( 'z\'' in filt):
            out.write(str(nbz) + ',' + 'z' ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbz) + ',' + 'z' ',' + '1' + '\n')

        if 'zcut' in filt:
            out.write(str(nbzcut) + ',' + 'zcut' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbzcut) + ',' + 'zcut' + ',' + '1' + '\n')

        if ('r' in filt) or( 'r\'' in filt):
            out.write(str(nbr) + ',' + 'r' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbr) + ',' + 'r' ',' + '1' + '\n')

        if ('i' in filt) or( 'i\'' in filt):
            out.write(str(nbi) + ',' + 'i' + ',' + '1' + '\n')
        else:
            out.write(str00 +  str(nbi) + ',' + 'i' ',' + '1' + '\n')

        if ('g' in filt) or( 'g\'' in filt):
            out.write(str(nbg) + ',' + 'g' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbg) + ',' + 'g' + ',' + '1' + '\n')

        if ('I+z' in filt) or( 'I+z\'' in filt):
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')
        else:
            out.write(str(nbIz) + ',' + 'I+z' + ',' + '1' + '\n')

        if 'Exo' in filt:
            out.write(str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbExo) + ',' + 'Exo' + ',' + '1' + '\n')

        if 'Clear' in filt:
            out.write(str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')
        else:
            out.write(str00 + str(nbClear) + ',' + 'Clear' + ',' + '1' + '\n')


def biasdark(t_now, Path, telescope, texps=None,bining_2=False):
    if telescope == 'Saint-Ex':
        with open(os.path.join(Path, str(t_now), 'Cal_biasdark'+ '_' +
                                                 datetime.strptime(t_now, '%Y-%m-%d').strftime('%m%d%Y')
                                                 + '.txt'), 'w') as out:
            str00 = ';'
            str1 = '#domeclose \n'
            str2 = '#nopreview \n'
            str3 = ' == bias dark exoplanet =='
            str4 = '#count '
            str5 = '#binning '
            str6 = '#interval '
            str7 = '#dark \n'
            str8 = '#shutdown \n'
            str9 = 'END'
            counts = ['9']*(len(texps)+1)
            binnings = ['1'] * (len(texps)+1)
            out.write(str1)
            out.write(str2)
            out.write(str00 + '\n')
            out.write(str00 + str3 + '\n')
            out.write(str00 + '\n')
            out.write(str4 + ','.join([str(counts) for counts in counts]) + '\n')
            out.write(str5 + ','.join([str(binnings) for binnings in binnings]) + '\n')
            out.write(str6 + '0,' + ','.join([str(texps) for texps in texps]) + '\n')
            out.write(str7)
            out.write(str8)
            out.write(str00 + '\n')
            out.write(str00 + str9 + '\n')
            out.write(str00 + '\n')
    elif telescope != 'Saint-Ex' and texps is not None:
        with open(os.path.join(Path, str(t_now), 'Cal_biasdark.txt'), 'w') as out:
            str00 = ';'
            str1 = '#domeclose \n'
            str2 = '#nopreview \n'
            str3 = ' == bias dark exoplanet =='
            str4 = '#count '
            str5 = '#binning '
            str6 = '#interval '
            str7 = '#dark \n'
            str8 = '#shutdown \n'
            str9 = 'END'
            counts = ['9']*(len(texps)+1)
            binnings = ['1'] * (len(texps)+1)
            out.write(str1)
            out.write(str2)
            out.write(str00 + '\n')
            out.write(str00 + str3 + '\n')
            out.write(str00 + '\n')
            out.write(str4 + ','.join([str(counts) for counts in counts]) + '\n')
            out.write(str5 + ','.join([str(binnings) for binnings in binnings]) + '\n')
            out.write(str6 + '0,' + ','.join([str(texps) for texps in texps]) + '\n')
            out.write(str7)
            out.write(str8)
            out.write(str00 + '\n')
            out.write(str00 + str9 + '\n')
            out.write(str00 + '\n')
    else:
        with open(os.path.join(Path, str(t_now), 'Cal_biasdark.txt'), 'w') as out:
            str00=';'
            str1 = '#domeclose \n'
            str2 = '#nopreview \n'
            str3 = ' == bias dark exoplanet =='
            str4 = '#count '
            str5 = '#binning '
            str6 = '#interval '
            str7 = '#dark \n'
            str8 = '#shutdown \n'
            str9 = 'END'
            if bining_2:
                texps = ['0', '15', '30', '60', '120', '0', '30', '60', '120', '240']
                counts = ['9'] * len(texps)
                binnings = ['1'] * (int(len(texps)/2)) + ['2'] * (int(len(texps)/2))
            else:
                texps = ['0', '15', '30', '60', '120']
                counts = ['9'] * len(texps)
                binnings = ['1'] * len(texps)
            out.write(str1)
            out.write(str2)
            out.write(str00 + '\n')
            out.write(str00 + str3 + '\n')
            out.write(str00 + '\n')
            out.write(str4 + ','.join([str(counts) for counts in counts]) + '\n')
            out.write(str5 + ','.join([str(binnings) for binnings in binnings]) + '\n')
            out.write(str6 + ','.join([str(texps) for texps in texps]) + '\n')
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


def haumea(t_now, date_start, date_end, count, filt, exptime,
           name_2, binning, Path, telescope,autofocus=True):
    waitlimit = 600
    date_start = np.datetime64(date_start)
    date_end = np.datetime64(date_end)
    hour0 = date_start.astype(object).hour
    hour0 = '{:02d}'.format(int(hour0))
    minute0 = date_start.astype(object).minute
    minute0 = '{:02d}'.format(int(minute0))
    hour1 = date_end.astype(object).hour
    hour1 = '{:02d}'.format(int(hour1))
    minute1 = date_end.astype(object).minute
    minute1 = '{:02d}'.format(int(minute1))
    with open(os.path.join(Path, str(t_now), 'Obj_haumea.txt'), 'w') as out:
        str00 = ';'
        str1 = '#waituntil 1, '
        str22 = '#domeopen\n'
        str2 = '#waitinlimits '
        str3 = '#nopreview \n'
        str4 = '#autofocus \n'
        str5 = '#count '
        str6 = '#binning '
        str7 = '#filter '
        str8 = '#interval '
        str9 = '#quitat '
        str10 = '#chain '
        for i in range(1, 2):
            out.write(str00 + '\n')
        out.write(str00 + ' ' + '136108 Haumea' + '\n')
        out.write(str00 + '\n')
        for i in range(1, 2):
            out.write(str00 + '\n')
        out.write(str1 + pd.to_datetime(str(date_start)).strftime(
            '%Y/%m/%d %H:%M:%S') + '\n')  # + str(hour0) + ':' + str(minute0) + '\n')
        out.write(str22)
        out.write('#chill -60\n')
        out.write(str2 + str(waitlimit) + '\n')
        out.write(str3)
        if not autofocus:
            out.write(str00 + str4)
        else:
            out.write(str4)
        out.write(str5 + str(count) + '\n')
        out.write(str6 + str(binning) + '\n')
        out.write(str7 + str(filt) + '\n')
        out.write(str8 + str(exptime) + '\n')
        out.write('D6108    0.26  0.15 K221L 218.88410  239.71462  122.12608   28.21248  0.1989312  0.00349648' +
                  '  42.9914375  1 MPO660785  3224  28 1955-2022 0.49 M-v 3Ek Pan        000A (136108) Haumea' + '\n')
        out.write(str(t_now) + '\n')
        out.write(str00 + '\n')
        out.write(str9 + str(hour1) + ':' + str(minute1) + '\n')
        if name_2 is None:
            out.write(str10 + 'Cal_flatdawn' + '.txt' + '\n')
        else:
            out.write(str10 + 'Obj_' + name_2 + '.txt' + '\n')
        out.write(str00 + '\n')


