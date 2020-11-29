import json
import orbits
import requests
import math
import mapbox
import numpy as np
import datetime as dt
import matplotlib.pyplot as plta
import julian as j
import pytz
from flask import Flask, render_template,request,jsonify, session
from flask_cors import CORS
import polyline as pl
import pandas as pd
coordinates = 0
coordinates2 = 0
utc = pytz.utc
RADIUS_EARTH = 6378
DEG_TO_RAD = lambda deg : deg * math.pi/180.0
from apscheduler.schedulers.background import BackgroundScheduler

geocoder = mapbox.Geocoder(access_token = 'pk.eyJ1IjoiYW5kcmV3Y3JhbXAiLCJhIjoiY2pvdmw4NzhoMThhczNrbzR4d2x0bGVhdyJ9.sc2tMk0EWnPkeCJWALbQ0g')
daysInMonth = [31, 29, 31, 30, 31, 30, 31, 31,  30, 31, 30, 31]
date = np.empty(50,dtype = 'object')
startTime = np.empty(50,dtype = 'object' )
startAzimuth = np.zeros(50)
endTime = np.zeros(50, dtype = 'object')
endAzimuth= np.zeros(50)
passes = 0
satPass = False

LATITUDECORRECTION = [["SE","SW"],
                        ["NE","NW"]]

def getCityCoordinates():
    coordinates = np.array([0.0,0.0])
    location = geocoder.forward('Welland,ON', limit = 1)
    location = location.json()
    coordinates[0] = location['features'][0]['center'][1]
    coordinates[1] = location['features'][0]['center'][0]
    return coordinates

def latLongToCartesian(latitude, longitude, radius):
    cartesian = np.array([0.0, 0.0, 0.0])
    cartesian[0] = radius*math.cos(latitude)*math.cos(longitude)
    cartesian[1] = radius*math.sin(longitude)*math.cos(latitude)
    cartesian[2] = radius*math.sin(latitude)
    return cartesian

def getDifferenceVector(observerVector, ISSVector):
    return np.subtract(ISSVector, observerVector)

def calculateElevationAngle(diffVector, earthVector):
    angle = math.acos(np.dot(earthVector,diffVector)/(np.linalg.norm(diffVector)*np.linalg.norm(earthVector)))

def lookAngle(city_coords, sat_coords):
    a = 6731.230
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
def getSemiMajorAxis(n):
    return (2.97554e15/((2*math.pi*n)**2))**(1.0/3.0)
     
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
        angle = lookAngle(observe_coord, sat_coord)
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
        check_passes(angle, i)
    print(count)
    lat.append(temp_lat)
    lon.append(temp_lon)
    temp_lat = []
    temp_lon = []
    for i in range(0, count+1):
        print(i)
        lonArray = np.array(lon[i])
        latArray = np.array(lat[i])
        coordinates.append(json.loads(pd.DataFrame(np.column_stack([lonArray, latArray])).to_json(orient='split'))['data']) 
    return coordinates

def check_passes(angle, index):
    global satPass
    global passes
    if  angle[0] > 3 and satPass != True:
        date_time = dt.datetime.now()
        date_time = date_time + dt.timedelta(hours=index*8.33e-4)
        date[passes] = date_time.strftime("%Y/%m/%d")
        startTime[passes] = date_time.strftime("%H:%M:%S")
        satPass = True
        startAzimuth[passes] = angle[1]
        startTimeTemp = float(index*8.33e-4)
    if satPass == True and angle[0] < 0:
        satPass = False
     #   passTime[passes] = i*8.33e-4-startTimeTemp
        endAzimuth[passes] = angle[1]
         #   if(passTime[passes] > 0.016666667): 
          #  passTime[passes] = passTime[passes] * 60
        date_time = dt.datetime.now()
        date_time = date_time + dt.timedelta(hours=index*8.33e-4)
        endTime[passes] = date_time.strftime("%H:%M:%S")
        passes = passes + 1
def getTimeFraction():
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

def getCoordIndex(start_time, time_step):
    current_time = j.to_jd(dt.datetime.now(utc), fmt='jd')
    index = int(math.floor((current_time - start_time) / time_step));
    return index

app = Flask(__name__)
app.secret_key = 'asdfasdfy4etyrejtryjfghnmmnb.,,.n,m.n,m.jkl.hj'
CORS(app, support_credentials=True)
@app.route("/", methods = ["GET","POST"])
def home():
    if session.get('satellite') == None:
        session['user_coord'] = [0.0,0.0]
        session['satellite_coord'] = [0.0,0.0]
        session['satellite'] = 25544
        session['tle'] = getTLE(session['satellite'])
        session['retrival_time'] = j.to_jd(dt.datetime.now(utc), fmt='jd')
    if request.method == "POST":
        new_coords = [float(request.form["latitude"]), float(request.form["longitude"])]
        new_id = request.form["satellite_id"]
        if new_id != session['satellite']  or new_coords != session['user_coord']:
            print(new_coords)
            session['satellite'] = new_id
            session['user_coord'] = new_coords
            session['tle'] = getTLE(session['satellite'])
            session['retrival_time'] = j.to_jd(dt.datetime.now(utc), fmt='jd')
    orbit_propogation = propogate_orbit(session['tle'], session['user_coord'])
    current_time = j.to_jd(dt.datetime.now(utc), fmt='jd')
    time_diff = current_time - session['retrival_time']
    session['satellite_coord'] = getFuturePosition(time_diff, session['tle'])
    session['look angle'] = lookAngle(session['user_coord'], session['satellite_coord'])
    print(orbit_propogation)
    return render_template("index.html",passes = passes,startTime = startTime, date = date,startAzimuth = startAzimuth,endTime = endTime,endAzimuth = endAzimuth, elev = round(session['look angle'][0],2), az = round(session['look angle'][1],2), lat = round(session['satellite_coord'][0],2), lon = round(session['satellite_coord'][1],2), sat_id = session['satellite'], latg = round(session['user_coord'][0],2), longg = round(session['user_coord'][1],2), coords = orbit_propogation, inclination = round(session['tle']['inclination']*180/math.pi,2), perigee = round(session['tle']['perigee']*180/math.pi,2), eccentricity = round(session['tle']['eccentricity'],5))


@app.route("/Update", methods = ["Get","POST"])
def update():
    lAngle = 0
    current_time = j.to_jd(dt.datetime.now(utc), fmt='jd')
    time_diff = current_time - session['retrival_time']
    session['satellite_coord'] = getFuturePosition(time_diff, session['tle'])
    lAngle = lookAngle(session['user_coord'],session['satellite_coord'])
    elevation = lAngle[0]
    azimuth = lAngle[1]
    return jsonify(elev = round(elevation,2), az = round(azimuth,2), lat = round(session['satellite_coord'][0],2), lon = round(session['satellite_coord'][1],2), latg = round(session['user_coord'][0],2), longg = round(session['user_coord'][1],2))

if (__name__ == "__main__"):
    app.run(host='localhost')
