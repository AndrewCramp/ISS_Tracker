import requests
import mapbox
geocoder = mapbox.Geocoder(access_token = 'pk.eyJ1IjoiYW5kcmV3Y3JhbXAiLCJhIjoiY2pvdmw4NzhoMThhczNrbzR4d2x0bGVhdyJ9.sc2tMk0EWnPkeCJWALbQ0g')

def getISSData():
    coordinates = []
    data = requests.get('http://api.open-notify.org/iss-now.json')
    if data.status_code == requests.codes.ok:
        print("Success\n")
    else:
        print("failure\n")
    data = data.json()
    coordinates[0] = data["iss_position"]["latitude"]
    coordinates[1] = data["iss_position"]["longitude"]
    return coordinates

def getCityCoordinates
    coordinates = []
    location = geocoder.forward('Kingston, ON', limit = 1)
    print(location.json()['features'][0]['center'])
    coordinates[0] = location.json()['features']['0']['center']
    coordinates[1] = location.json()['features']['0']['center']
    return coordinates

getISSData()
