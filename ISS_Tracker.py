import json
import requests
import math
import mapbox
import numpy as np
import datetime as dt
import matplotlib.pyplot as plta
import julian as j
import pytz
from flask import Flask, render_template,request,jsonify
import polyline as pl
import pandas as pd
coordinates = 0
coordinates2 = 0
utc = pytz.utc
RADIUS_EARTH = 6378
ISS_ORBIT = 408+RADIUS_EARTH
semiMajorAxis = 6796000
geocoder = mapbox.Geocoder(access_token = 'pk.eyJ1IjoiYW5kcmV3Y3JhbXAiLCJhIjoiY2pvdmw4NzhoMThhczNrbzR4d2x0bGVhdyJ9.sc2tMk0EWnPkeCJWALbQ0g')
inclination = 0.0
ascendingNode = 0.0
perigee = 0.0
epochTime = 0.0
epochAnomaly = 0.0
meanMotion = 0.0
daysInMonth = [31, 29, 31, 30, 31, 30, 31, 31,  30, 31, 30, 31]
observeLat = 0
observeLon = 0


LATITUDECORRECTION = [["SE","SW"],
                        ["NE","NW"]]
def getISSData():
    coordinates = np.array([0.0,0.0])
    data = requests.get('http://api.open-notify.org/iss-now.json')
    data = data.json()
    coordinates[0] = data["iss_position"]["latitude"]
    coordinates[1] = data["iss_position"]["longitude"]
    return coordinates

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

def lookAngle(cityLat, cityLong, satLat, satLong):
    a = 6731.230
    angle = np.array([0.0,0.0])
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


def getTLE():
    global epochTime
    global epochAnomaly
    global inclination
    global perigee
    global ascendingNode
    global meanMotion
    global eccentricity
    data = requests.get('https://www.n2yo.com/rest/v1/satellite/tle/25544&apiKey=C3LD4X-SJY6XX-HDLDBW-4DIZ')
    tle = [' ',' ']
    temp = data.json()['tle'].split('\r\n')
    tle[0] = temp[0]
    tle[1] = temp[1]
    inclination = float(tle[1][9:16])*math.pi/180.0
    ascendingNode = float(tle[1][17:25])*math.pi/180.0
    perigee = float(tle[1][34:42])*math.pi/180.0
    epochTime = float(tle[0][20:32])
    epochAnomaly = float(tle[1][43:51])*math.pi/180.0
    meanMotion = float(tle[1][52:63])
    eccentricity = float(tle[1][26:33])/10000000

def getMeanAnomaly(deltat):
    n = meanMotion
    M0 = epochAnomaly
    M =  M0 + 2*math.pi*(n*deltat-int(n*deltat)-int((M0+2*math.pi*(n*deltat-int(n*deltat)))/(2*math.pi)))
    return M
def getEccentricAnomaly(M):
    iterations =5
    error = 0
    e = eccentricity
    i = 0 
    E = M if(e < 0.8) else math.pi
    F = E - e * math.sin(M)-M
    while((abs(F) > error) and (i < iterations)):
        E = E -F/(1.0-e*math.cos(E))
        F - E-e*math.sin(E) - M
        i+=1
    return E


def getTrueAnomaly(E):
    e = eccentricity
    ta =  math.acos((math.cos(E)-e)/(1-e*math.cos(E)))
    ta2 = 2*math.pi-ta
    if(abs(ta2-E) > abs(ta-E)):
        return ta
    else:
        return ta2
def getSemiMajorAxis():
    n = meanMotion
    return (2.97554e15/((2*math.pi*n)**2))**(1.0/3.0)
     
def getPerigeeDistance(a):
    e = eccentricity
    return a*(1-e)

def getRadius(perigee, trueAnomaly):
    e = eccentricity
    p = perigee
    v = trueAnomaly
    return (p*(1+e))/(1+e*math.cos(v))

def getRAPrecession(a, deltat):
    Re = 1.0
    a1 = a/RADIUS_EARTH
    i = inclination
    n = meanMotion
    e = eccentricity
    J2 = 1.0826e-3
    d1temp = 3*J2*Re**2*(3*math.cos(i)**2-1)
    d1temp2 = 4*(a1**2)*(1-e)**(3.0/2.0)
    d1 = (d1temp/d1temp2)
    a0 = -a1*(134*d1**3/81 + d1**2 +d1/3 -1)
    p0 = a0*(1-e**2)
    r = ascendingNode
    return r + 2*math.pi*(-3*J2*Re**2*n*deltat*math.cos(inclination)/(2*p0**2))

def getPerigeePrecession(a, deltat):
    Re = 1.0
    a1 = a/RADIUS_EARTH
    i = inclination
    n = meanMotion
    e = eccentricity
    J2 = 1.0826e-3
    d1temp = 3*J2*Re**2*(3*math.cos(i)**2-1)
    d1temp2 = 4*(a1**2)*(1-e)**(3.0/2.0)
    d1 = (d1temp/d1temp2)
    a0 = -a1*(134*d1**3/81 + d1**2 +d1/3 -1)
    p0 = a0*(1-e**2)
    omega = perigee
    return omega + 2*math.pi*(3*J2*Re**2*n*deltat*(5*math.cos(inclination)**2-1)/(4*p0**2))

def getArgumentofLatitude(perigeeP, v):
    return perigeeP + v - 2*math.pi*(int((perigeeP+v)/(2*math.pi)))

def getRADifference(mu):
    i = inclination
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

def getFuturePosition(hours):
    days = hours/24
    deltat = getTimeFraction()-epochTime+days
    meanAnomaly = getMeanAnomaly(deltat)
    E = getEccentricAnomaly(meanAnomaly)
    v = getTrueAnomaly(E)
    a = getSemiMajorAxis()
    p = getPerigeeDistance(a)
    RAP = getRAPrecession(a, deltat)
    perigeeP = getPerigeePrecession(a, deltat)
    mu = getArgumentofLatitude(perigeeP, v)
    RAdiff = getRADifference(mu)
    r = getRadius(p,v)
    alphag = getGeocentricRA(RAdiff, RAP)
    deltag = getGeocentricDeclination(mu,RAdiff)
    X1 = a*(math.cos(E)-eccentricity)
    Y1 = a*(math.sqrt(1-eccentricity**2)*math.sin(E))
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

def plotLine():
    global coordinates
    global coordinates2
    lat = []
    lon = []
    lat2 = []
    lon2 = []
    switch = 0
    tempp = 5000
    for i in range(1800):
        temp = getFuturePosition(i*8.33e-4)
        if abs(temp[1] - tempp) > 15 and tempp != 5000:
            switch = 1
        if switch == 0:
            lat.append(float(temp[0]))
            lon.append(float(temp[1]))
            tempp = temp[1]
        else:
            lat2.append(float(temp[0]))
            lon2.append(float(temp[1]))
    lonArray = np.array(lon)
    latArray = np.array(lat)
    lon2Array = np.array(lon2)
    lat2Array = np.array(lat2)
    line = np.column_stack([lonArray, latArray])
    line2 = np.column_stack([lon2Array, lat2Array])
    coordinates = (json.loads(pd.DataFrame(line).to_json(orient='split'))['data']) 
    coordinates2 = (json.loads(pd.DataFrame(line2).to_json(orient='split'))['data']) 
    
    

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

def plotLat():
    j = 0
    count = 0
    xAxis = np.empty([63])
    yAxis = np.empty([63])
    yAxis2 = np.empty([63])
    while(j<2*math.pi):
        meanAnomaly = j
        xAxis[count] = meanAnomaly
        eccentricAnomaly = 0.0
        semiMajorAxis = 6796000
        iterations =500
        error = 0
        i = 0 
        E = meanAnomaly if(eccentricity < 0.8) else math.pi
        F = E - eccentricity * math.sin(meanAnomaly)-meanAnomaly
        while((abs(F) > error) and (i < iterations)):
            E = E -F/(1.0-eccentricity*math.cos(E))
            F - E-eccentricity*math.sin(E) - meanAnomaly
            i+=1
        X1 = semiMajorAxis*(math.cos(E)-eccentricity)
        Y1 = semiMajorAxis*(math.sqrt(1-eccentricity**2)*math.sin(E))
        x = math.cos(perigee) * X1 - math.sin(perigee) * Y1
        y = math.sin(perigee) * X1 + math.cos(perigee) * Y1
        z = math.sin(inclination) * x
        x = math.cos(inclination) * x
        xtemp = x
        x = math.cos(ascendingNode) * xtemp - math.sin(ascendingNode) * y
        y = math.sin(ascendingNode) * xtemp +math.cos(ascendingNode) * y
        r= math.sqrt(x**2+y**2+z**2)
        lat = math.asin(z/semiMajorAxis)
        yAxis[count] = lat*180.0/math.pi

    
        coord1 = np.array([[X1],[Y1],[0.0]])
        
        AzPerigee = np.array([[math.cos(perigee),math.sin(perigee),0.0],
                          [-1*math.sin(perigee),math.cos(perigee),0.0],
                          [0.0,0.0,1]])
        AzAscending = np.array([[math.cos(ascendingNode),math.sin(ascendingNode),0.0],
                          [-1*math.sin(ascendingNode),math.cos(ascendingNode),0.0],
                          [0.0,0.0,1]])


        Ax = np.array([[1,0,0],
                   [0,math.cos(inclination),math.sin(inclination)],
                   [0,-1*math.sin(inclination),math.cos(inclination)]])

        temp1 = AzAscending.transpose().dot(Ax.transpose()).dot(AzPerigee.transpose())
        cartesian = temp1.dot(coord1)
        r= math.sqrt(cartesian[0]**2+cartesian[1]**2+cartesian[2]**2)
        lat = math.atan2(cartesian[2],(math.sqrt(cartesian[0]**2+cartesian[1]**2)))
        yAxis2[count] = lat*180/math.pi
        j+=0.1
        count = count + 1
    plt.figure()
    plt.plot(xAxis,yAxis)
    plt.figure()
    plt.plot(xAxis, yAxis2)
    plt.show()

app = Flask(__name__)
@app.route("/", methods = ["GET","POST"])
def home():
    global observLat
    global observeLon
    elevation = 0
    azimuth = 0
    lAngle = 0
    ISSCoord = np.array([0.0,0.0])
    cityCoord = [0.0,0.0]
    ISSCoord = getISSData()
    ISSCoord = ISSCoord*(math.pi/180)
    getTLE()
    location = getFuturePosition(0)
    latitude = location[0]
    longitude = location[1]
    plotLine()
    if request.method == "POST":
        cityCoord = [float(request.form["latitude"])*math.pi/180, float(request.form["longitude"])*math.pi/180]
        observeLat = float(request.form["latitude"])*math.pi/180
        observeLon = float(request.form["longitude"])*math.pi/180
        lAngle = lookAngle(cityCoord[0],cityCoord[1],ISSCoord[0],ISSCoord[1])
        elevation = lAngle[0]
        azimuth = lAngle[1]
    return render_template("index.html", elev = round(elevation,5), az = round(azimuth,5), lat = round(latitude,5), lon = round(longitude,5), latg = round(cityCoord[0]*180/math.pi,5), longg = round(cityCoord[1]*180/math.pi,5), coords = coordinates,coords2 = coordinates2, inclination = round(inclination*180/math.pi,5), perigee = round(perigee*180/math.pi,5), eccentricity = round(eccentricity,5))


@app.route("/Update", methods = ["Get","POST"])
def update():
    elevation = 0
    azimuth = 0
    lAngle = 0
    ISSCoord = np.array([0.0,0.0])
    cityCoord = [0.0,0.0]
    ISSCoord = getISSData()
    ISSCoord = ISSCoord*(math.pi/180)
    latitude = ISSCoord[0]*180/math.pi
    longitude = ISSCoord[1]*180/math.pi
    cityCoord = [observeLat,observeLon]
    lAngle = lookAngle(cityCoord[0],cityCoord[1],ISSCoord[0],ISSCoord[1])
    elevation = lAngle[0]
    azimuth = lAngle[1]
    print(lAngle)
    return jsonify(elev = round(elevation,5), az = round(azimuth,5), lat = round(latitude,5), lon = round(longitude,5), latg = round(cityCoord[0]*180/math.pi,5), longg = round(cityCoord[1]*180/math.pi,5))

if (__name__ == "__main__"):
    app.run()
