import json
from orbits import orbits
import math
import mapbox
import numpy as np
import datetime as dt
import julian as j
from flask import Flask, render_template,request,jsonify, session
from flask_cors import CORS
from timezonefinder import TimezoneFinder
import pytz
tf = TimezoneFinder()
geocoder = mapbox.Geocoder(access_token = 'pk.eyJ1IjoiYW5kcmV3Y3JhbXAiLCJhIjoiY2pvdmw4NzhoMThhczNrbzR4d2x0bGVhdyJ9.sc2tMk0EWnPkeCJWALbQ0g')

def getCityCoordinates():
    coordinates = np.array([0.0,0.0])
    location = geocoder.forward('Welland,ON', limit = 1)
    location = location.json()
    coordinates[0] = location['features'][0]['center'][1]
    coordinates[1] = location['features'][0]['center'][0]
    return coordinates



app = Flask(__name__)
app.secret_key = 'asdfasdfy4etyrejtryjfghnmmnb.,,.n,m.n,m.jkl.hj'
CORS(app, support_credentials=True)
@app.route("/", methods = ["GET","POST"])
def home():
    retrival_time_local = dt.datetime.now(orbits.utc)
    if session.get('satellite') == None:
        session['user_coord'] = [0.0,0.0]
        session['satellite_coord'] = [0.0,0.0]
        session['satellite'] = 25544
        session['tle'] = orbits.getTLE(session['satellite'])
        timezone_str = tf.timezone_at(lng = session['user_coord'][1],lat = session['user_coord'][0])
        timezone = pytz.timezone(timezone_str)
        session['retrival_time'] = j.to_jd(dt.datetime.now(orbits.utc), fmt='jd')
        retrival_time_local= dt.datetime.now(timezone)
    if request.method == "POST":
        new_coords = [float(request.form["latitude"]), float(request.form["longitude"])]
        new_id = request.form["satellite_id"]
        if new_id != session['satellite']  or new_coords != session['user_coord']:
            session['satellite'] = new_id
            session['user_coord'] = new_coords
            timezone_str = tf.timezone_at(lng = session['user_coord'][1],lat = session['user_coord'][0])
            timezone = pytz.timezone(timezone_str)
            session['tle'] = orbits.getTLE(session['satellite'])
            session['retrival_time'] = j.to_jd(dt.datetime.now(orbits.utc), fmt='jd')
            retrival_time_local= dt.datetime.now(timezone)
    orbit_propogation = orbits.propogate_orbit(session['tle'], session['user_coord'])
    temp_coords = []
    for coord_array in orbit_propogation:
        temp_coords.extend(coord_array)
    current_time = j.to_jd(dt.datetime.now(orbits.utc), fmt='jd')
    time_diff = current_time - session['retrival_time']
    sat_passes = orbits.check_passes(temp_coords, session['user_coord'], orbits.getSemiMajorAxis(session['tle']['mean motion']), retrival_time_local)
    startTime = []
    startAzimuth = []
    endTime = []
    endAzimuth = []
    endTime = []
    date = []
    startTime = sat_passes['start time']
    startAzimuth = sat_passes['start azimuth']
    endTime = sat_passes['end time']
    endAzimuth = sat_passes['end azimuth']
    date = sat_passes['date']
    passes = sat_passes['passes']
    session['satellite_coord'] = orbits.getFuturePosition(time_diff, session['tle'])
    session['look angle'] = orbits.lookAngle(session['user_coord'], session['satellite_coord'],orbits.getSemiMajorAxis(session['tle']['mean motion']))
    return render_template("index.html",passes = passes,startTime = startTime, date = date,startAzimuth = startAzimuth,endTime = endTime,endAzimuth = endAzimuth, elev = round(session['look angle'][0],2), az = round(session['look angle'][1],2), lat = round(session['satellite_coord'][0],2), lon = round(session['satellite_coord'][1],2), sat_id = session['satellite'], latg = round(session['user_coord'][0],2), longg = round(session['user_coord'][1],2), coords = orbit_propogation, inclination = round(session['tle']['inclination']*180/math.pi,2), perigee = round(session['tle']['perigee']*180/math.pi,2), eccentricity = round(session['tle']['eccentricity'],5))


@app.route("/Update", methods = ["Get","POST"])
def update():
    lAngle = 0
    current_time = j.to_jd(dt.datetime.now(orbits.utc), fmt='jd')
    time_diff = current_time - session['retrival_time']
    session['satellite_coord'] = orbits.getFuturePosition(time_diff, session['tle'])
    lAngle = orbits.lookAngle(session['user_coord'],session['satellite_coord'], orbits.getSemiMajorAxis(session['tle']['mean motion']))
    elevation = lAngle[0]
    azimuth = lAngle[1]
    return jsonify(elev = round(elevation,2), az = round(azimuth,2), lat = round(session['satellite_coord'][0],2), lon = round(session['satellite_coord'][1],2), latg = round(session['user_coord'][0],2), longg = round(session['user_coord'][1],2))

if (__name__ == "__main__"):
    app.run(host='localhost')
