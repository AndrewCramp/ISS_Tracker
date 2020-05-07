import requests
import math
import mapbox
import numpy as np
import datetime as dt
import matplotlib.pyplot as plta
import julian as j
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


LATITUDECORRECTION = [["SE","SW"],
                        ["NE","NW"]]
def getISSData():
    coordinates = np.array([0.0,0.0])
    data = requests.get('http://api.open-notify.org/iss-now.json')
    if data.status_code == requests.codes.ok:
        print("Success\n")
    else:
        print("failure\n")
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
    print(angle)

def lookAngle(cityLat, cityLong, satLat, satLong):
    a = getSemiMajorAxis()
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
    print(sector)
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
    print(temp)
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
    print("M: ",M*180.0/math.pi)
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
    return 2*math.pi-ta
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
    print(perigee)
    return omega + 2*math.pi*(3*J2*Re**2*n*deltat*(5*math.cos(inclination)**2-1)/(4*p0**2))

def getArgumentofLatitude(perigeeP, v):
    return perigeeP + v - 2*math.pi*(int((perigeeP+v)/(2*math.pi)))

def getRADifference(mu):
    i = inclination
    if(0<i and i < math.pi/2 and 0 < mu and mu <math.pi):
        return math.acos(math.cos(mu)/(1-math.sin(i)**2*math.sin(mu)**2)**0.5)
    else:
        return 2*math.pi - math.acos(math.cos(mu)/(1-math.sin(i)**2*math.sin(mu)**2)**0.5)
def getGeocentricRA(RAdiff, asscendingNodeP):
    omega = asscendingNodeP
    return RAdiff + omega - 2*math.pi*(int(RAdiff+omega/(2*math.pi)))

def getGeocentricDeclination(mu,RAdiff):
    x = -1 if(math.sin(mu) < 0) else 1
    print(math.sin(mu))
    print(x)
    return x*math.acos(math.cos(mu)/math.cos(RAdiff))

def getFuturePosition():
    deltat = getTimeFraction()-epochTime+0.00347222222
    print(deltat)
    print("time: ",getTimeFraction())
    print("Mean anomaly at epoch: ",epochAnomaly*180/math.pi)
    meanAnomaly = getMeanAnomaly(deltat)
    print("mean anomaly: ", meanAnomaly*180/math.pi)
    semiMajorAxis = 6796000
    E = getEccentricAnomaly(meanAnomaly)
    print("Eccentric Anomaly: ",E*180/math.pi)
    print("True anomaly: ",getTrueAnomaly(E)*180/math.pi)
    v = getTrueAnomaly(E)
    print("semi: ", getSemiMajorAxis())
    a = getSemiMajorAxis()
    print("perigee distance: ", getPerigeeDistance(a))
    p = getPerigeeDistance(a)
    print("radius: ", getRadius(p,v))
    print("RA: ", getRAPrecession(a, deltat)*180/math.pi)
    RAP = getRAPrecession(a, deltat)
    print("pergee: ", getPerigeePrecession(a, deltat)*180/math.pi)
    perigeeP = getPerigeePrecession(a, deltat)
    mu = getArgumentofLatitude(perigeeP, v)
    print("argument of Latitude: ", mu*180/math.pi)
    print("ra diff: ", getRADifference(mu)*180/math.pi)
    RAdiff = getRADifference(mu)
    r = getRadius(p,v)
    alphag = getGeocentricRA(RAdiff, RAP)
    deltag = getGeocentricDeclination(mu,RAdiff)
    print("GeocentricDeclination: ", deltag)
    X1 = a*(math.cos(E)-eccentricity)
    Y1 = a*(math.sqrt(1-eccentricity**2)*math.sin(E))
   # x = math.cos(perigee) * X1 - math.sin(perigee) * Y1
   # y = math.sin(perigee) * X1 + math.cos(perigee) * Y1
   # z = math.sin(inclination) * x
   # x = math.cos(inclination) * x
   # xtemp = x
   # x = math.cos(ascendingNode) * xtemp - math.sin(ascendingNode) * y
   # y = math.sin(ascendingNode) * xtemp +math.cos(ascendingNode) * y
    
    x = r*math.cos(alphag)*math.cos(deltag)
    y = r*math.sin(alphag)*math.cos(deltag)
    z = r*math.sin(deltag)

    latg = 42.99*math.pi/180
    longg = -79.24*math.pi/180
    siderealt = 211.25*math.pi/180
    rg =((math.cos(latg)/RADIUS_EARTH)**2+(math.sin(latg)/6356.72)**2)**(-1.0/2.0)
    print(rg)
    xg = rg*math.cos(siderealt)*math.cos(latg)
    yg = rg*math.sin(siderealt)*math.cos(latg)
    zg = rg*math.sin(latg)
    xs = x-xg
    ys = y-yg
    zs = z-zg
    rs = math.sqrt(xs**2+ys**2+zs**2)
    print(x)
    print(y)
    print(z)
    r=math.sqrt(x**2+y**2+z**2)
    print(r)
    lat = math.asin(z/r)
    lon = math.atan2(y, x)
    if(lon < 0):
        lon += 2*math.pi
    print("LAT:")
    print(lat*180/math.pi)
    print(lon*180/math.pi)


    coord1 = np.array([[X1],[Y1],[0.0]])
    AzPerigee = np.array([[math.cos(perigee),math.sin(perigee),0.0],
                          [-1.0*math.sin(perigee),math.cos(perigee),0.0],
                          [0.0,0.0,1.0]])
    AzAscending = np.array([[math.cos(ascendingNode),math.sin(ascendingNode),0.0],
                          [-1.0*math.sin(ascendingNode),math.cos(ascendingNode),0.0],
                          [0.0,0.0,1.0]])


    Ax = np.array([[1.0,0.0,0.0],
                   [0.0,math.cos(inclination),math.sin(inclination)],
                   [0.0,-1.0*math.sin(inclination),math.cos(inclination)]])
    temp1 = AzAscending.transpose().dot(Ax.transpose()).dot(AzPerigee.transpose())
    cartesian = temp1.dot(coord1)
    r= math.sqrt(cartesian[0]**2+cartesian[1]**2+cartesian[2]**2)
    v = math.acos((math.cos(E)-eccentricity)/(1-eccentricity*math.cos(E)))
    lat = math.asin(cartesian[2]/r)
    lon = math.atan2(cartesian[1], cartesian[0])
    if(lon < 0):
        lon += 2*math.pi
    print("LAT:")
    print(lat*180/math.pi)
    print(lon*180/math.pi)

def getTimeFraction():
    currentDT = dt.datetime.now()
    month = currentDT.month
    days = 0 
    for i in range(month-1):
        days = days + daysInMonth[i]
    days = days + currentDT.day
    seconds = (currentDT.hour+4)*60*60
    seconds = seconds + currentDT.minute*60
    seconds = seconds + currentDT.second
    fraction = seconds/86400
    return days+fraction

def getLMST():
    jd = j.to_jd(dt.datetime.now())-2451545
    T = jd/36525.0
    GMST = 24110.54841+8640184.812866*T+0.093104*T**2-0.0000062*T**3
    GMST = GMST*0.000277778
    LMST = GMST + 280.76
    print(GMST)

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

        
ISSCoord = np.array([0.0,0.0])
cityCoord = np.array([0.0,0.0])
ISSCoord = getISSData()
cityCoord = getCityCoordinates()
ISSCoord = ISSCoord*(math.pi/180)
cityCoord = cityCoord*(math.pi/180)
getTLE()
print(lookAngle(cityCoord[0],cityCoord[1],ISSCoord[0],ISSCoord[1]))
getFuturePosition()
getLMST()
