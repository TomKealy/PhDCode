import json 
import urllib

response = urllib.urlopen("http://search.twitter.com/search.json?q=microsoft")

js = json.load(response)
