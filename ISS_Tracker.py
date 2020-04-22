import requests

def getISSData():
    data = requests.get('http://api.open-notify.org/iss-now.json')
    if data.status_code == requests.codes.ok:
        print("Success")
    else:
        print("failure")

getISSData()
