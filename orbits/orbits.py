import json
import settings
import requests
import math
import numpy as np
import datetime as dt
import pandas as pd
import pytz
import numpy as np
import julian as j
utc = pytz.utc
RADIUS_EARTH = 6378
DEG_TO_RAD = lambda deg : deg * math.pi/180.0

def lookAngle(city_coords, sat_coords, a):
    LATITUDECORRECTION = [["SE","SW"],["NE","NW"]]
    cityLat = DEG_TO_RAD(city_coords[0])
    cityLong = DEG_TO_RAD(city_coords[1])
    satLat =  DEG_TO_RAD(sat_coords[0])
    satLong = DEG_TO_RAD(sat_coords[1])
    angle = [0.0,0.0]
    gamma = math.acos(math.sin(satLat)*math.sin(cityLat)+math.cos(satLat)*math.cos(cityLat)*math.cos(satLong-cityLong))
    elevation = math.acos(math.sin(gamma)/(math.sqrt(1+math.pow(RADIUS_EARTH/a,2)-2.000*(RADIUS_EARTH/a)*math.cos(gamma))))
    elevation = elevation * 180.000/math.pi
    if(gamma > math.acos(RADIUS_EARTH/a)):
        elevation = elevation * -1
    alpha = math.asin(math.sin(math.fabs(cityLong-satLong))*math.cos(satLat)/math.sin(gamma))
    north = 0
    west = 0
    if(satLong < cityLong):
        west = 1
    x = (math.tan(cityLat)*math.cos(satLong-cityLong))
    y = math.tan(satLat)
    if(x < y):
        north = 1
    azimuth = getAzimuth(LATITUDECORRECTION[north][west],alpha*180/math.pi)
    angle[0] = elevation
    angle[1] = azimuth
    return angle


def getAzimuth(sector, alpha):
    if(sector == "NW"):
        return 360-alpha
    elif(sector == "NE"):
        return alpha
    elif(sector == "SW"):
        return 180 + alpha
    elif(sector == "SE"):
        return 180 - alpha


def getTLE(satellite_id):
    data = requests.get(f'https://data.ivanstanojevic.me/api/tle/{satellite_id}')
    tle = [' ',' ']
    tle[0] = data.json()['line1']
    tle[1] = data.json()['line2']
    inclination = float(tle[1][9:16])*math.pi/180.0
    ascendingNode = float(tle[1][17:25])*math.pi/180.0
    perigee = float(tle[1][34:42])*math.pi/180.0
    epochTime = float(tle[0][19:32])
    epochAnomaly = float(tle[1][43:51])*math.pi/180.0
    meanMotion = float(tle[1][52:63])
    eccentricity = float(tle[1][26:33])/10000000
    tle_data = {
            'inclination': inclination,
            'ascending node': ascendingNode,
            'perigee': perigee,
            'epoch time': epochTime,
            'epoch anomaly': epochAnomaly,
            'mean motion': meanMotion,
            'eccentricity': eccentricity
        }
    return tle_data

def getMeanAnomaly(deltat,M0, n):
    M =  M0 + 2*math.pi*(n*deltat-int(n*deltat)-int((M0+2*math.pi*(n*deltat-int(n*deltat)))/(2*math.pi)))
    return M
def getEccentricAnomaly(M,e):
    iterations =5
    error = 0
    i = 0 
    E = M if(e < 0.8) else math.pi
    F = E - e * math.sin(M)-M
    while((abs(F) > error) and (i < iterations)):
        E = E -F/(1.0-e*math.cos(E))
        F - E-e*math.sin(E) - M
        i+=1
    return E


def getTrueAnomaly(E, e):
    ta =  math.acos((math.cos(E)-e)/(1-e*math.cos(E)))
    ta2 = 2*math.pi-ta
    if(abs(ta2-E) > abs(ta-E)):
        return ta
    else:
        return ta2

def getPerigeeDistance(a, e):
    return a*(1-e)

def getRadius(perigee, trueAnomaly, e, p):
    v = trueAnomaly
    return (p*(1+e))/(1+e*math.cos(v))

def getRAPrecession(a, deltat, i, n, e, r):
    Re = 1.0
    a1 = a/RADIUS_EARTH
    J2 = 1.0826e-3
    d1temp = 3*J2*Re**2*(3*math.cos(i)**2-1)
    d1temp2 = 4*(a1**2)*(1-e)**(3.0/2.0)
    d1 = (d1temp/d1temp2)
    a0 = -a1*(134*d1**3/81 + d1**2 +d1/3 -1)
    p0 = a0*(1-e**2)
    return r + 2*math.pi*(-3*J2*Re**2*n*deltat*math.cos(i)/(2*p0**2))

def getPerigeePrecession(a, deltat, i, n, e, perigee):
    Re = 1.0
    a1 = a/RADIUS_EARTH
    J2 = 1.0826e-3
    d1temp = 3*J2*Re**2*(3*math.cos(i)**2-1)
    d1temp2 = 4*(a1**2)*(1-e)**(3.0/2.0)
    d1 = (d1temp/d1temp2)
    a0 = -a1*(134*d1**3/81 + d1**2 +d1/3 -1)
    p0 = a0*(1-e**2)
    omega = perigee
    return omega + 2*math.pi*(3*J2*Re**2*n*deltat*(5*math.cos(i)**2-1)/(4*p0**2))

def getArgumentofLatitude(perigeeP, v):
    return perigeeP + v - 2*math.pi*(int((perigeeP+v)/(2*math.pi)))

def getRADifference(mu, i):
    if(0<i and i < math.pi/2 and 0 < mu and mu <math.pi or math.pi/2 < i and i < math.pi and math.pi < mu and mu < math.pi*2):
        return math.acos(math.cos(mu)/(1-math.sin(i)**2*math.sin(mu)**2)**0.5)
    else:
        return 2*math.pi - math.acos(math.cos(mu)/(1-math.sin(i)**2*math.sin(mu)**2)**0.5)

def getGeocentricRA(RAdiff, asscendingNodeP):
    omega = asscendingNodeP
    return RAdiff + omega - 2*math.pi*(int(RAdiff+omega/(2*math.pi)))

def getGeocentricDeclination(mu,RAdiff):
    x = -1 if(math.sin(mu) < 0) else 1
    return x*math.acos(math.cos(mu)/math.cos(RAdiff))

def getSemiMajorAxis(n):
    G = 6.67408e-11 
    M = 5.972e24
    a = ((G*M / (2*math.pi*n**2)) ** (1.0/3.0))
    return a



def getFuturePosition(hours, tle):
    days = hours/24
    deltat = getTimeFraction()-tle['epoch time']+days
    meanAnomaly = getMeanAnomaly(deltat, tle['epoch anomaly'],tle['mean motion'])
    E = getEccentricAnomaly(meanAnomaly, tle['eccentricity'])
    v = getTrueAnomaly(E, tle['eccentricity'])
    a = getSemiMajorAxis(tle['mean motion'])
    p = getPerigeeDistance(a,tle['eccentricity'])
    RAP = getRAPrecession(a, deltat, tle['inclination'], tle['mean motion'], tle['eccentricity'], tle['ascending node'])
    perigeeP = getPerigeePrecession(a, deltat, tle['inclination'], tle['mean motion'], tle['eccentricity'], tle['perigee'])
    mu = getArgumentofLatitude(perigeeP, v)
    RAdiff = getRADifference(mu, tle['inclination'])
    r = getRadius(p,v, tle['eccentricity'], tle['perigee'])
    alphag = getGeocentricRA(RAdiff, RAP)
    deltag = getGeocentricDeclination(mu,RAdiff)
    X1 = a*(math.cos(E)-tle['eccentricity'])
    Y1 = a*(math.sqrt(1-tle['eccentricity']**2)*math.sin(E))
    x = r*math.cos(alphag)*math.cos(deltag)
    y = r*math.sin(alphag)*math.cos(deltag)
    z = r*math.sin(deltag)
    r=math.sqrt(x**2+y**2+z**2)
    lat = math.asin(z/r)
    lon = math.atan2(y, x)*180.00/math.pi
    adjustment = getGMST(hours)*15
    lon = lon - adjustment
    while(lon < 0):
        lon = lon + 360
    pos = [lat*180/math.pi, lon]
    return pos

def propogate_orbit(tle, observe_coord):
    lat = []
    lon = []
    lat2 = []
    lon2 = []
    coordinates = []
    temp_lat = []
    temp_lon = []
    count = 0
    tempp = 5000
    for i in range(36000):
        sat_coord = getFuturePosition(i*8.33e-4, tle)
        sat_coord.reverse()
        angle = lookAngle(observe_coord, sat_coord, getSemiMajorAxis(tle['mean motion']))
        sat_coord.reverse()
        if i < 36000:
            if abs(sat_coord[1] - tempp) > 15 and tempp != 5000:
                count = count + 1
                lat.append(temp_lat)
                lon.append(temp_lon)
                temp_lat = []
                temp_lon = []
                tempp = sat_coord[1]
            else:
                temp_lat.append(float(sat_coord[0]))
                temp_lon.append(float(sat_coord[1]))
                tempp = sat_coord[1]
    lat.append(temp_lat)
    lon.append(temp_lon)
    temp_lat = []
    temp_lon = []
    for i in range(0, count+1):
        lonArray = np.array(lon[i])
        latArray = np.array(lat[i])
        coordinates.append(json.loads(pd.DataFrame(np.column_stack([lonArray, latArray])).to_json(orient='split'))['data']) 
    return coordinates

def check_passes(coordinates, city_coords, a, retrival_time):
    satPass = False
    passes = 0
    pass_list = []
    count = -1
    sat_pass = {
        'date': [],
        'start time': [],
        'end time': [],
        'start azimuth': [],
        'end azimuth': [],
    }
    for coord in coordinates:
        count = count + 1
        coord.reverse()
        angle = lookAngle(city_coords, coord, a)
        coord.reverse()
        if  angle[0] > 0 and satPass != True:
            date_time = retrival_time
            date_time = date_time + dt.timedelta(hours=count*8.33e-4)
            sat_pass['date'].append(date_time.strftime("%Y/%m/%d"))
            sat_pass['start time'].append(date_time.strftime("%H:%M:%S"))
            sat_pass['start azimuth'].append(angle[1])
            satPass = True
        if satPass == True and angle[0] < 0:
            satPass = False
            sat_pass['end azimuth'].append(angle[1])
            date_time = retrival_time
            date_time = date_time + dt.timedelta(hours=count*8.33e-4)
            sat_pass['end time'].append(date_time.strftime("%H:%M:%S"))
            passes = passes + 1
    sat_pass['passes'] = passes
    return sat_pass

def getTimeFraction():
    daysInMonth = [31, 29, 31, 30, 31, 30, 31, 31,  30, 31, 30, 31]
    currentDT = dt.datetime.now(utc)
    month = currentDT.month
    days = 0 
    for i in range(month-1):
        days = days + daysInMonth[i]
    days = days + currentDT.day
    seconds = (currentDT.hour)*60*60
    seconds = seconds + currentDT.minute*60
    seconds = seconds + currentDT.second
    fraction = seconds/86400
    return days+fraction

def getGMST(deltat):
    jd = j.to_jd(dt.datetime.now(utc), fmt='jd')
    midnight = math.floor(jd)+0.5
    daysSinceMidnight = jd - midnight
    hoursSinceMidnight = daysSinceMidnight*24
    daysSinceEpoch = jd - 2451545
    centuriesSinceEpoch = daysSinceEpoch / 36525
    wholeDaysSinceEpoch = midnight - 2451545.0
    GMST = (6.697374558
    + 0.06570982441908 * wholeDaysSinceEpoch
    + 1.00273790935 * hoursSinceMidnight
    + 0.000026 * centuriesSinceEpoch**2)
    hours = GMST%24
    return hours+deltat
