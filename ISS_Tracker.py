import requests
import math
import mapbox
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
RADIUS_EARTH = 63781000
ISS_ORBIT = 408000+RADIUS_EARTH
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
    angle = np.array([0.0,0.0])
    gamma = math.acos(math.sin(satLat)*math.sin(cityLat)+math.cos(satLat)*math.cos(cityLat)*math.cos(satLong-cityLong))
    elevation = math.acos(math.sin(gamma)/(math.sqrt(1+math.pow(RADIUS_EARTH/ISS_ORBIT,2)-2.000*(RADIUS_EARTH/ISS_ORBIT)*math.cos(gamma))))
    elevation = elevation * 180.000/math.pi
    if(gamma > math.acos(RADIUS_EARTH/ISS_ORBIT)):
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
    print(north)
    print(west)
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
    print("anomaly: ") 
    print(epochAnomaly)
    print(epochTime)
    meanMotion = float(tle[1][52:63])*2*math.pi
    eccentricity = float(tle[1][26:33])/10000000
    print(eccentricity)


def getFuturePosition():
    print(getTimeFraction()-epochTime)
    meanAnomaly = epochAnomaly + meanMotion*(getTimeFraction()+0.00347-epochTime)
    while(meanAnomaly>2*math.pi):
        meanAnomaly -= 2*math.pi
        print(meanAnomaly)
    print("motion:")
    print(meanMotion)
    print(meanAnomaly)
    eccentricAnomaly = 0.0
    semiMajorAxis = 6796000
    iterations =5
    error = 0
    i = 0 
    E = meanAnomaly if(eccentricity < 0.8) else math.pi
    print(E)
    F = E - eccentricity * math.sin(meanAnomaly)-meanAnomaly
    while((abs(F) > error) and (i < iterations)):
            E = E -F/(1.0-eccentricity*math.cos(E))
            F - E-eccentricity*math.sin(E) - meanAnomaly
            i+=1
            print(E)
    X1 = semiMajorAxis*(math.cos(E)-eccentricity)
    Y1 = semiMajorAxis*(math.sqrt(1-eccentricity**2)*math.sin(E))
    x = math.cos(perigee) * X1 - math.sin(perigee) * Y1
    y = math.sin(perigee) * X1 + math.cos(perigee) * Y1
    z = math.sin(inclination) * x
    x = math.cos(inclination) * x
    xtemp = x
    x = math.cos(ascendingNode) * xtemp - math.sin(ascendingNode) * y
    y = math.sin(ascendingNode) * xtemp +math.cos(ascendingNode) * y
    print(x)
    print(y)
    print(z)

    r= math.sqrt(x**2+y**2+z**2)
    print(r)
    lat = math.asin(z/r)
    lon = math.atan2(y, x)
    if(lon < 0):
        lon += 2*math.pi
    print(lat*180/math.pi)
    print(lon*180/math.pi)


    coord1 = np.array([[X1],[Y1],[0.0]])
    print(coord1)
    AzPerigee = np.array([[math.cos(perigee),math.sin(perigee),0.0],[0.0,math.cos(perigee),math.sin(perigee)],[0.0,-1*math.sin(perigee),math.cos(perigee)]])
    AzAscending = np.array([[math.cos(ascendingNode),math.sin(ascendingNode),0.0],[0.0,math.cos(ascendingNode),math.sin(ascendingNode)],[0.0,-1*math.sin(ascendingNode),math.cos(ascendingNode)]])
    Ax = np.array([[1,0,0],[0,math.cos(inclination),math.sin(inclination)],[0,-1*math.sin(inclination),math.cos(inclination)]])
    print(Ax)
    print(AzPerigee)
    print(AzAscending)
    temp1 = AzAscending.transpose().dot(Ax.transpose()).dot(AzPerigee.transpose())
    print(temp1)
    cartesian = temp1.dot(coord1)
    print(cartesian)
    r= math.sqrt(cartesian[0]**2+cartesian[1]**2+cartesian[2]**2)
    print(r)
    lat = math.asin(cartesian[2]/r)
    lon = math.atan2(cartesian[1], cartesian[0])
    if(lon < 0):
        lon += 2*math.pi
    print(lat*180/math.pi)
    print(lon*180/math.pi)

def getTimeFraction():
    currentDT = dt.datetime.now()
    month = currentDT.month
    days = 0 
    for i in range(month-1):
        days = days + daysInMonth[i]
    days = days + currentDT.day
    seconds = currentDT.hour*60*60
    seconds = seconds + currentDT.minute*60
    seconds = seconds + currentDT.second
    fraction = seconds/86400
    return days+fraction


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
        iterations =5
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
        AzPerigee = np.array([[math.cos(perigee),math.sin(perigee),0.0],[0.0,math.cos(perigee),math.sin(perigee)],[0.0,-1*math.sin(perigee),math.cos(perigee)]])
        AzAscending = np.array([[math.cos(ascendingNode),math.sin(ascendingNode),0.0],[0.0,math.cos(ascendingNode),math.sin(ascendingNode)],[0.0,-1*math.sin(ascendingNode),math.cos(ascendingNode)]])
        Ax = np.array([[1,0,0],[0,math.cos(inclination),math.sin(inclination)],[0,-1*math.sin(inclination),math.cos(inclination)]])
        temp1 = AzAscending.transpose().dot(Ax.transpose()).dot(AzPerigee.transpose())
        cartesian = temp1.dot(coord1)
        r= math.sqrt(cartesian[0]**2+cartesian[1]**2+cartesian[2]**2)
        lat = math.asin(cartesian[2]/r)
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
print(lookAngle(cityCoord[0],cityCoord[1],ISSCoord[0],ISSCoord[1]))
getTLE()
getFuturePosition()
plotLat()
