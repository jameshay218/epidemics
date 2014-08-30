# Import the modules
import requests
import json

# Make it a bit prettier..
print "-" * 30
print "This will show the Most Popular Videos on YouTube"
print "-" * 30

# Get the feed
r = requests.get("https://www.googleapis.com/youtube/analytics/v1/reports?ids=channel%3D%3DUCO1M1qbs3Y2Xx-vaydzFPuQ&start-date=2014-05-01&end-date=2014-06-30&metrics=views&key={YOUR_API_KEY}")
r.text
# Convert it to a Python dictionary
data = json.loads(r.text)

# Loop through the result. 
for item in data['data']['items']:
    print item
    
