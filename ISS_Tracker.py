import requests
import math
import mapbox
import numpy as np
RADIUS_EARTH = 63781000000
ISS_ORBIT = 408000000+RADIUS_EARTH
geocoder = mapbox.Geocoder(access_token = 'pk.eyJ1IjoiYW5kcmV3Y3JhbXAiLCJhIjoiY2pvdmw4NzhoMThhczNrbzR4d2x0bGVhdyJ9.sc2tMk0EWnPkeCJWALbQ0g')

LATITUDECORRECTION = [["NW","NE"],
                        ["SE","SW"]]
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
    location = geocoder.forward('Kingston,ON', limit = 1)
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
    elevation = math.acos(math.sin(gamma)/(math.sqrt(1+math.pow(RADIUS_EARTH/ISS_ORBIT,2)-2*(RADIUS_EARTH/ISS_ORBIT)*math.cos(gamma))))
    elevation = elevation * 180/math.pi
    if(gamma > math.acos(RADIUS_EARTH/ISS_ORBIT)):
        elevation = elevation * -1
    alpha = math.asin(math.sin(math.fabs(cityLong-satLong))*math.cos(satLat)/math.sin(gamma))
    north = 0
    west = 0
    if(satLong < cityLong):
        west = 1
    x = (math.tan(cityLat)*math.cos(satLong-cityLong))
    y = math.tan(satLat)
    if(x < 7):
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


ISSCoord = np.array([0.0,0.0])
cityCoord = np.array([0.0,0.0])
ISSCoord = getISSData()
cityCoord = getCityCoordinates()
ISSCoord = ISSCoord*(math.pi/180)
cityCoord = cityCoord*(math.pi/180)
print(lookAngle(cityCoord[0],cityCoord[1],ISSCoord[0],ISSCoord[1]))
